/*
 * File: ImageTemplate.cs
 *
 * Author: Eric Peterson
 * Eric.J.Peterson@gmail.com
 * University of California, Riverside
 * Smart Tools Lab 2009
 * 
 * Use at your own risk.  This code is not maintained and not guaranteed to work.
 * We take no responsibility for any harm this code may cause.
 */

using System;
using System.Collections.Generic;
using System.Text;
using System.Diagnostics;
using Sketch;
using Utilities;

namespace ImageAligner
{
    /// <summary>
    /// Template for a symbol which is compared using an Image-Based Recognizer.
    /// The Image-Based Recognizer compares individual strokes to between two
    /// templates. These individual strokes are sized and placed by a bounding 
    /// box around the entire group of strokes being compared. The ink in one 
    /// unknown stroke can only be compared to one known template stroke.
    /// 
    /// Strokes are mapped between an unknown template and a known template. When
    /// the Recognize function is called, all likely mappings between the templates
    /// are made and then evaluated. In order to perform and evaluate the mapping 
    /// some pre-processing must be done by a single-stroke recognizer (SSR), such 
    /// as Rubine's or the Dollar Recognizer, to determine the mostly likely 
    /// component a stroke belongs to (e.g. a stroke could be recognized as being
    /// the front arc of an AND gate, or the shaft of an arrow).
    /// 
    /// When recognizing, it will return an ImageTemplateResult, which includes 
    /// references to the unknown and best-matching template, as well as the 
    /// score of this match. The score ranges from -infinity (very bad) to 1.0 (very 
    /// good), however in practical use it ranges from -0.2 to 1.0 (the only reason 
    /// it can go below 0.0 is if there are errors found in the match). The result 
    /// also contains information about any errors that are encountered in mapping 
    /// the templates, for instance if there is a bubble stroke in the template that 
    /// is not found in the unknown, the ImageTemplateResult will report that there is 
    /// a Missing Bubble. Another example would be if there were an extra stroke in the 
    /// unknown (perhaps the stroke is actually a wire instead of a gate - but looks 
    /// like a TopArc), in which case it could report Extra TopArc.
    /// 
    /// Example: the unknown shape is classified as a gate, and contains two strokes:
    /// a BackArc and a TopArc (parts of an OR gate). For the sake of this example, 
    /// let's assume that these two strokes were correctly recognized by the SSR as 
    /// being a BackArc and a TopArc. This unknown shape is made into a template and
    /// compared against a list of pre-defined templates, including a 3-stroke OR gate 
    /// comprised of a BackArc, TopArc, and BottomArc. The most likely mapping between 
    /// the unknown and the template is BackArc:BackArc, TopArc:TopArc, and 
    /// UnMapped:BottomArc. The image-based recognizer will then compare the ink in for 
    /// each of these pairs, reporting an image score composed of Haussdorff and Modified 
    /// Haussdorf distances. Let's say that the BackArc:BackArc comparison had a Haussdorf
    /// distance of 3.0 and a Modified Haussdorf distance of 1.7, and the TopArc:TopArc
    /// mapping had distances of 2.0 and 1.3. The combined result (ImageTemplateResult) 
    /// would have a Haussdorff distance of 3.0 (max of all individual symbol distances)
    /// and a Modified Haussdorf distance of 1.5 (average of all individual symbol distances).
    /// The Image score for this mapping would be calculated by:
    /// 1.0 / (1.0 + Haussdorf * ScaleFactor1 + ModifiedHaussdorf * ScaleFactor2)
    /// where the scale factors are based on the size of the bitmaps to give more meaningful 
    /// score ranges (e.g. a score of 0.9 is achievable and good using the scale factors, whereas
    /// without them a good score might be 0.15 while the max is still 1.0). Currently the scale
    /// factors are 1/48 for the Haussdorf distance and 2/48 for the Modified Haussdorf distance.
    /// This Image only score ignores the fact that there is an unmapped stroke - the overall 
    /// score penalizes for each unmapped stroke (currently set to 0.1 / stroke). 
    /// Thus the image only score would be 0.89, and the overall score would be 0.79.
    /// </summary>
    [DebuggerDisplay("{m_SymbolInfo.SymbolType}: {m_ShapeParts.Count} Strokes (Definition Template: {m_IsFullyKnown})")]
    [Serializable]
    public class ImageTemplate
    {
        #region Constants and Parameters

        /// <summary>
        /// File containing a pre-determined confusion matrix
        /// </summary>
        string c_ConfusionMatrixFile = "Code\\Recognition\\ImageAligner\\allConfusion.txt";

        /// <summary>
        /// This is the minimum allowable probability for a symbol's mapping,
        /// in other words, if we are trying to map a shape part X to
        /// another shape part Y, and the probability is less than this
        /// number, that mapping will not be considered.
        /// </summary>
        const float MIN_PROBABILITY = 0.01f;

        /// <summary>
        /// Once we have scored each of the mappings we can sort them based on their
        /// score. This number tells us how many of the top mappings to try and 
        /// recognize. Essentially this should be a trade-off between speed and accuracy.
        /// </summary>
        const int NUM_TOP_MAPPINGS = 5;
        
        /// <summary>
        /// This parameter works with the NUM_TOP_MAPPINGS to determine
        /// whether a mapping should be considered. If the mapping score
        /// is greater than this number, it means that the mapping is 
        /// pretty bad, and has a low chance of being a good match for the 
        /// template, so we shouldn't try to recognize with that mapping.
        /// </summary>
        const float MAX_ALLOWABLE_MAPPING_SCORE = 3.0f;

        /// <summary>
        /// Fudge factor when determining whether a stroke's bounding box
        /// is to the left/right/above/below of another stroke's bbox
        /// </summary>
        const float MAX_BBOX_INCLUSION = 0.3f;

        #endregion


        #region Member Variables

        /// <summary>
        /// All the general information about the shape.
        /// User info, shape descriptions, under what conditions the symbol was drawn.
        /// </summary>
        SymbolInfo m_SymbolInfo;

        /// <summary>
        /// Information from the Single-Stroke Recognizer
        /// </summary>
        Dictionary<Substroke, string> m_ShapeParts;

        /// <summary>
        /// All ImageSymbols, key is a hash of substroke IDs and the bounding box's hash code
        /// </summary>
        Dictionary<int, ImageSymbol> m_Symbols;

        /// <summary>
        /// Stores the bounding boxes for various combinations of Substrokes
        /// int key is a hash of the substroke guid hashcodes with ^
        /// </summary>
        Dictionary<int, System.Drawing.Rectangle> m_Boxes;

        /// <summary>
        /// The complete shape that this template was created from
        /// </summary>
        Shape m_CompleteShape;

        /// <summary>
        /// Confusion matrix values for the single stroke recognizer (SSR)
        /// </summary>
        ConfusionMatrix m_SSRConfusionMatrix;

        /// <summary>
        /// Essentially tells whether the template is known to be correct
        /// </summary>
        bool m_IsFullyKnown = false;

        Dictionary<Substroke, Dictionary<string, List<Substroke>>> m_StrokeLocationRelations;

        #endregion


        #region Constructors

        /// <summary>
        /// Creates a new ImageTemplate for an unknown shape
        /// The strokes are labeled with beliefs from a single-stroke 
        /// recognizer in the shapeParts variable
        /// </summary>
        /// <param name="completeShape">Shape object containing all of the strokes</param>
        /// <param name="shapeParts">Believed parts of the shape from the single-stroke recognizer</param>
        public ImageTemplate(Shape completeShape, Dictionary<Substroke, string> shapeParts)
        {
            InitilizeUnknown(completeShape, shapeParts);

            GetConfusionMatrix(c_ConfusionMatrixFile);
        }

        /// <summary>
        /// Creates a new ImageTemplate from a shape that has been correctly 
        /// labeled. These templates are to be used for recognition later
        /// </summary>
        /// <param name="completeShape">Shape object containing all of the strokes</param>
        /// <param name="shapeParts">Known parts of the shape</param>
        /// <param name="info">General information about the symbol, such as user or platform</param>
        public ImageTemplate(Shape completeShape, Dictionary<Substroke, string> shapeParts, SymbolInfo info)
        {
            InitializeKnown(completeShape, shapeParts, info);

            GetConfusionMatrix(c_ConfusionMatrixFile);
        }

        /// <summary>
        /// Creates a new ImageTemplate for an unknown shape
        /// The strokes are labeled with beliefs from a single-stroke 
        /// recognizer in the shapeParts variable
        /// </summary>
        /// <param name="completeShape">Shape object containing all of the strokes</param>
        /// <param name="shapeParts">Believed parts of the shape from the single-stroke recognizer</param>
        /// <param name="matrix">Confusion matrix for SSR</param>
        public ImageTemplate(Shape completeShape, Dictionary<Substroke, string> shapeParts, ConfusionMatrix matrix)
        {
            InitilizeUnknown(completeShape, shapeParts);

            m_SSRConfusionMatrix = matrix;
        }

        /// <summary>
        /// Creates a new ImageTemplate from a shape that has been correctly 
        /// labeled. These templates are to be used for recognition later
        /// </summary>
        /// <param name="completeShape">Shape object containing all of the strokes</param>
        /// <param name="shapeParts">Known parts of the shape</param>
        /// /// <param name="info">General information about the symbol, such as user or platform</param>
        /// <param name="matrix">Confusion matrix for SSR</param>
        public ImageTemplate(Shape completeShape, Dictionary<Substroke, string> shapeParts, SymbolInfo info, ConfusionMatrix matrix)
        {
            InitializeKnown(completeShape, shapeParts, info);

            m_SSRConfusionMatrix = matrix;
        }

        /// <summary>
        /// Constructor helper to avoid duplicate code
        /// </summary>
        /// <param name="completeShape">Shape object containing all of the strokes</param>
        /// <param name="shapeParts">Known parts of the shape</param>
        /// <param name="info">General information about the symbol, such as user or platform</param>
        private void InitializeKnown(Shape completeShape, Dictionary<Substroke, string> shapeParts, SymbolInfo info)
        {
            m_SymbolInfo = info;
            m_ShapeParts = shapeParts;
            m_CompleteShape = completeShape;
            m_Symbols = new Dictionary<int, ImageSymbol>();
            m_Boxes = new Dictionary<int, System.Drawing.Rectangle>();
            m_IsFullyKnown = true;

            List<Substroke> strokes = completeShape.SubstrokesL;
            System.Drawing.Rectangle bbox = Compute.BoundingBox(strokes.ToArray());
            m_Boxes.Add(GetHashStrokes(strokes), bbox);
            int hashShape = GetHashStrokes(strokes) ^ GetHashBbox(bbox);
            m_Symbols.Add(hashShape, new ImageSymbol(strokes, bbox, m_SymbolInfo));

            // Create an ImageSymbol for each of the strokes, use the bounding box for 
            // the entire shape
            foreach (Substroke s in strokes)
            {
                ImageSymbol symbol = new ImageSymbol(s.PointsAsSysPoints, bbox, info);
                int hash = s.Id.GetHashCode() ^ GetHashBbox(bbox);
                m_Symbols.Add(hash, symbol);
            }

            m_StrokeLocationRelations = GetLocationRelations(completeShape);
        }

        /// <summary>
        /// Creates a new ImageTemplate for an unknown shape
        /// The strokes are labeled with beliefs from a single-stroke 
        /// recognizer in the shapeParts variable
        /// </summary>
        /// <param name="completeShape">Shape object containing all of the strokes</param>
        /// <param name="shapeParts">Believed parts of the shape from the single-stroke recognizer</param>
        private void InitilizeUnknown(Shape completeShape, Dictionary<Substroke, string> shapeParts)
        {
            m_SymbolInfo = new SymbolInfo();
            m_Boxes = GetGroupBoxes(GetStrokeCombinations(completeShape.SubstrokesL));
            m_ShapeParts = shapeParts;
            m_Symbols = new Dictionary<int, ImageSymbol>();
            m_CompleteShape = completeShape;

            m_StrokeLocationRelations = GetLocationRelations(completeShape);
        }

        /// <summary>
        /// Initializes the SSR confusion matrix with a made-up matrix
        /// </summary>
        private void InitializeConfusionMatrix()
        {
            List<string> labels = General.leafLabels;
            
            m_SSRConfusionMatrix = new ConfusionMatrix(labels);

            foreach (string label in labels)
                for (int i = 0; i < 50; i++)
                    m_SSRConfusionMatrix.Add(label, label);

            m_SSRConfusionMatrix.IncreaseAllByOne();
        }

        /// <summary>
        /// Gets a confusion matrix from file
        /// </summary>
        /// <param name="filename"></param>
        private void GetConfusionMatrix(string filename)
        {
            string filePath = System.IO.Directory.GetCurrentDirectory();
            if (filePath.Contains("\\Code\\"))
                filePath = filePath.Substring(0, filePath.IndexOf("\\Code\\") + 1);
            else if (filePath.Contains("\\Sketch\\"))
                filePath = filePath.Substring(0, filePath.IndexOf("\\Sketch\\") + 8);
            else if (filePath.Contains("\\Trunk\\"))
                filePath = filePath.Substring(0, filePath.IndexOf("\\Trunk\\") + 7);
            filePath += c_ConfusionMatrixFile;
            m_SSRConfusionMatrix = new ConfusionMatrix(General.leafLabels);
            m_SSRConfusionMatrix.LoadFromFile(filePath);
        }

        #endregion


        #region Recognition

        public ImageTemplateResult Recognize(List<ImageTemplate> templates)
        {
            List<ImageTemplateResult> results = Recognize(templates, 1);

            if (results.Count == 0)
                return null;
            else
                return results[0];
        }

        /// <summary>
        /// Recognizes the current unknown template using a list of known templates.
        /// Determines the ImageTemplateResult for each template in the list and
        /// returns the best match. If no viable match is found, the function will
        /// return null.
        /// </summary>
        /// <param name="templates"></param>
        /// <returns></returns>
        public List<ImageTemplateResult> Recognize(List<ImageTemplate> templates, int n)
        {
            List<ImageTemplateResult> results = new List<ImageTemplateResult>();

            foreach (ImageTemplate template in templates)
            {
                ImageTemplateResult result = Recognize(template);

                if (result != null)
                    results.Add(result);
            }

            results.Sort(delegate(ImageTemplateResult R1, ImageTemplateResult R2) { return R1.Score.CompareTo(R2.Score); });
            results.Reverse();

            List<ImageTemplateResult> output = new List<ImageTemplateResult>(n);
            for (int i = 0; i < n; i++)
                if (results.Count > i)
                    output.Add(results[i]);

            return output;
        }

        /// <summary>
        /// Gets the recognition result (ImageTemplateResult) between this (unknown) 
        /// and the template (known). This is accomplished by attempting to find the
        /// best mapping between the two templates. This mapping is used to determine 
        /// the image score where ink in the unknown's stroke A can only be compared 
        /// to the ink in the template's stroke 1 (if the mapping is A:1).
        /// </summary>
        /// <param name="template"></param>
        /// <returns></returns>
        public ImageTemplateResult Recognize(ImageTemplate template)
        {
            if (template.m_StrokeLocationRelations == null)
                template.m_StrokeLocationRelations = template.GetLocationRelations(template.m_CompleteShape);

            // Find the most probable matchings: this.ImageSymbols --> template.ImageSymbols
            SortedList<double, StrokeMapping> bestMappings = GetMappingsTree(template);

            ImageTemplateResult bestResult = null;
            int n = 1;
            foreach (KeyValuePair<double, StrokeMapping> kv in bestMappings)
            {
                double mappingScore = kv.Key;
                StrokeMapping map = kv.Value;

                if (n > NUM_TOP_MAPPINGS || mappingScore >= MAX_ALLOWABLE_MAPPING_SCORE)
                    break;
                n++;

                // Gets the actual image recognizer score between the two templates using 
                // the current mapping between strokes
                ImageTemplateResult result = GetImageScore(map, this, template);

                if (result != null)
                    if (bestResult == null || (result.Score != double.NaN && result.Score > bestResult.Score))
                        bestResult = result;
            }

            return bestResult;
        }

        /// <summary>
        /// This function goes through each of the stroke mappings between the templates
        /// and computes the Haussdorf and Modified Haussdorf distances between the symbols.
        /// This function then creates an ImageTemplateResult which computes the composite 
        /// score based on the individual results.
        /// </summary>
        /// <param name="map">Stroke mapping indicating which stroke each of the unknown 
        /// strokes get compared to</param>
        /// <param name="uShape">Image Template for the unknown shape</param>
        /// <param name="tShape">Image Template for the known template shape</param>
        /// <returns>Result of matching the unknown to the template using this mapping</returns>
        private ImageTemplateResult GetImageScore(StrokeMapping map, ImageTemplate uShape, ImageTemplate tShape)
        {
            Dictionary<ImageSymbol, ImageSymbolResult> results = new Dictionary<ImageSymbol, ImageSymbolResult>();
            System.Drawing.Rectangle ubbox, tbbox;

            // Get the bounding box for all mapped strokes in the unknown.
            // This bounding box is used scale and position the individual strokes.
            List<Substroke> uGroup = map.StrokeKeys;
            int hash1 = GetHashStrokes(uGroup);
            if (uShape.m_Boxes.ContainsKey(hash1))
                ubbox = uShape.m_Boxes[hash1];
            else
                ubbox = Utilities.Compute.BoundingBox(uGroup.ToArray());

            // Get the bounding box for all mapped strokes in the known template.
            // If all strokes have been mapped, this box has already been computed.
            List<Substroke> tGroup = new List<Substroke>(tShape.m_CompleteShape.SubstrokesL);
            int hash2 = GetHashStrokes(tGroup);
            if (tShape.m_Boxes.ContainsKey(hash2))
                tbbox = tShape.m_Boxes[hash2];
            else
                tbbox = Utilities.Compute.BoundingBox(tGroup.ToArray());

            // Go through each stroke-stroke mapping and find image result
            foreach (KeyValuePair<Substroke, Substroke> kv in map.Pairs)
            {
                ImageSymbol unknown, template;

                // Get the ImageSymbol for the unknown, create it if necessary
                Substroke uStroke = kv.Key;
                int uHash = uStroke.Id.GetHashCode() ^ GetHashBbox(ubbox);
                if (uShape.m_Symbols.ContainsKey(uHash))
                    unknown = uShape.m_Symbols[uHash];
                else
                {
                    unknown = new ImageSymbol(uStroke.PointsAsSysPoints, ubbox, uShape.m_SymbolInfo);
                    uShape.m_Symbols.Add(uHash, unknown);
                }

                // Get the ImageSymbol for the template, create it if necessary
                Substroke tStroke = kv.Value;
                int tHash = tStroke.Id.GetHashCode() ^ GetHashBbox(tbbox);
                if (tShape.m_Symbols.ContainsKey(tHash))
                    template = tShape.m_Symbols[tHash];
                else
                {
                    template = new ImageSymbol(tStroke.PointsAsSysPoints, tbbox, tShape.m_SymbolInfo);
                    tShape.m_Symbols.Add(tHash, template);
                }

                // Recognize the unknown
                ImageSymbolResult result = unknown.Recognize(template);
                results.Add(unknown, result);
            }

            // Create the composite result using the individual ImageSymbolResults
            ImageTemplateResult resultTemplate = new ImageTemplateResult(uShape, tShape, map, results);

            return resultTemplate;
        }

        /// <summary>
        /// Searches for the best mappings between the unknown and template strokes
        /// </summary>
        /// <param name="template">template to get mappings with</param>
        /// <returns>Sorted results based on mapping score (best first)</returns>
        private SortedList<double, StrokeMapping> GetMappingsTree(ImageTemplate template)
        {
            List<StrokeMapping> mappings = new List<StrokeMapping>();
            SortedList<double, StrokeMapping> sortedMappings = new SortedList<double, StrokeMapping>();

            // Recurse and get all possible combinations
            GetMappingsTreeRecurse(this.Strokes, template.Strokes, new StrokeMapping(), ref mappings, 
                this, template, this.m_SSRConfusionMatrix);
            
            // Score each map and insert into the sorted list
            foreach (StrokeMapping mapping in mappings)
            {
                if (mapping.Count == 0) continue;
                
                double score = GetMappingScore(mapping, this, template);
                
                if (score == 0.0)
                    score += 0.000000001;

                while (sortedMappings.ContainsKey(score))
                    score += score / 1000000.0;

                sortedMappings.Add(score, mapping);
            }

            return sortedMappings;
        }

        /// <summary>
        /// Recursive function which builds stroke-mappings between 2 templates
        /// </summary>
        /// <param name="unknowns">Remaining strokes in the unknown template</param>
        /// <param name="templates">Remaining strokes in the known template</param>
        /// <param name="map">Current mapping of strokes</param>
        /// <param name="all">List of all mappings made so far</param>
        /// <param name="uParts">Dictionary of unknown substrokes to their believed type</param>
        /// <param name="tParts">Dictionary of template substrokes to their known type</param>
        /// <param name="matrix">Confusion matrix to use for rejecting very bad mappings</param>
        private void GetMappingsTreeRecurse(List<Substroke> unknowns, List<Substroke> templates, 
            StrokeMapping map, ref List<StrokeMapping> all, 
            ImageTemplate unknown, ImageTemplate template,
            ConfusionMatrix matrix)
        {
            // Add a valid map
            if (map.Count > 0 && !all.Contains(map))
                all.Add(map);

            // End recursion if there are no remaining strokes in
            // the unknown or the template
            if (unknowns.Count == 0 || templates.Count == 0)
                return;

            Dictionary<Substroke, string> uParts = unknown.m_ShapeParts;
            Dictionary<Substroke, string> tParts = template.m_ShapeParts;

            // Branch the search tree from each unknown stroke
            foreach (Substroke u in unknowns)
            {
                // Get the remaining unknown strokes - remove the current stroke
                List<Substroke> rU = new List<Substroke>(unknowns);
                rU.Remove(u);

                // Branch the search tree from each template stroke
                foreach (Substroke t in templates)
                {
                    // Get the remaining template strokes - remove the current stroke
                    List<Substroke> rT = new List<Substroke>(templates);
                    rT.Remove(t);

                    // Get the probability of this particular stroke-stroke mapping pair
                    double prob = matrix.GetProbability(uParts[u], tParts[t]);

                    // Ignore this mapping pair if it is lower than a given threshold
                    if (prob < MIN_PROBABILITY)
                        continue;

                    // Determine whether this match makes sense spatially
                    foreach (KeyValuePair<Substroke, Substroke> kvp in map.Pairs)
                    {
                        Substroke tStroke = kvp.Value;
                        Substroke uStroke = kvp.Key;
                        if (template.m_StrokeLocationRelations[t]["Left"].Contains(tStroke)
                            && !unknown.m_StrokeLocationRelations[u]["Left"].Contains(uStroke))
                            continue;

                        if (template.m_StrokeLocationRelations[t]["Right"].Contains(tStroke)
                            && !unknown.m_StrokeLocationRelations[u]["Right"].Contains(uStroke))
                            continue;

                        if (template.m_StrokeLocationRelations[t]["Above"].Contains(tStroke)
                            && !unknown.m_StrokeLocationRelations[u]["Above"].Contains(uStroke))
                            continue;

                        if (template.m_StrokeLocationRelations[t]["Below"].Contains(tStroke)
                            && !unknown.m_StrokeLocationRelations[u]["Below"].Contains(uStroke))
                            continue;
                    }

                    // Clone the current map so that we have a new object
                    StrokeMapping newMap = new StrokeMapping(map);

                    // Add the current mapping pair to the newly cloned map
                    newMap.Add(u, t);

                    // Continue the recursion with the new remaining strokes and new map
                    GetMappingsTreeRecurse(rU, rT, newMap, ref all, unknown, template, matrix);
                }
            }
        }

        /// <summary>
        /// Score the current mapping, penalize for unmapped strokes.
        /// Lower number --> Better mapping
        /// Score is the average of the absolute value of the log 
        /// probability of each mapping pair.
        /// </summary>
        /// <param name="mapping">Particular mapping to score</param>
        /// <param name="unknown">Unknown template</param>
        /// <param name="template">Known template</param>
        /// <returns>Mapping Score</returns>
        private double GetMappingScore(StrokeMapping mapping, ImageTemplate unknown, ImageTemplate template)
        {
            Dictionary<Substroke, string> uParts = unknown.m_ShapeParts;
            Dictionary<Substroke, string> tParts = template.m_ShapeParts;

            double score = 0.0;

            foreach (KeyValuePair<Substroke, Substroke> kv in mapping.Pairs)
            {
                if (!uParts.ContainsKey(kv.Key) || !tParts.ContainsKey(kv.Value))
                    throw new Exception("Substroke is unknown in the template");

                double prob = m_SSRConfusionMatrix.GetProbability(uParts[kv.Key], tParts[kv.Value]);

                // Should never need this, but just in case - make sure 
                // the probability is at least the MIN_PROBABILITY (especially non-zero)
                prob = Math.Max(prob, MIN_PROBABILITY);

                // Add the log probability to the total score - this is averaged later
                score += Math.Abs(Math.Log(prob));
            }

            // Will be used for averaging later
            int count = mapping.Count;

            // The probability for unmapped strokes is the same as the minimum
            // probability that can be obtained between any stroke mapping
            double unmappedPenalty = MAX_ALLOWABLE_MAPPING_SCORE;

            // Go through each stroke in the unknowns, penalizing for each
            // unmapped stroke
            foreach (Substroke s in unknown.Strokes)
                if (!mapping.ContainsKey(s))
                {
                    score += unmappedPenalty;
                    count++;
                }

            // Go through each stroke in the tempates, penalizing for each
            // unmapped stroke
            foreach (Substroke s in template.Strokes)
                if (!mapping.ContainsValue(s))
                {
                    score += unmappedPenalty;
                    count++;
                }

            // Fills things out in the mapping object
            mapping.DetermineUnmapped(unknown, template);

            // Make sure we won't be dividing by zero
            if (count == 0)
                return 1000.0;

            return score / count;
        }

        #endregion


        #region Helpers

        /// <summary>
        /// Gets all combinations of strokes for a given list
        /// </summary>
        /// <param name="strokes"></param>
        /// <returns></returns>
        private Dictionary<int, List<Substroke>> GetStrokeCombinations(List<Substroke> strokes)
        {
            Dictionary<int, List<Substroke>> allCombinations = new Dictionary<int, List<Substroke>>();
            List<Substroke> copy = new List<Substroke>(strokes);

            // Recurse to find all combinations
            GetStrokeCombinations(copy, ref allCombinations);

            return allCombinations;
        }

        /// <summary>
        /// Recursive function which searches for all combinations of remaining substrokes
        /// </summary>
        /// <param name="strokes">All remaining substrokes</param>
        /// <param name="all">List of all combinations so far</param>
        private void GetStrokeCombinations(List<Substroke> strokes, ref Dictionary<int, List<Substroke>> all)
        {
            if (strokes.Count == 0)
                return;

            int hash = GetHashStrokes(strokes);

            if (!all.ContainsKey(hash))
                all.Add(hash, strokes);

            foreach (Substroke s in strokes)
            {
                // Clone the list of substrokes
                List<Substroke> copy = new List<Substroke>(strokes);
                
                // Remove the current stroke
                copy.Remove(s);

                // Recurse over remaing strokes
                GetStrokeCombinations(copy, ref all);
            }
        }

        /// <summary>
        /// Computes bounding boxes for a set of substroke lists.
        /// Bounding boxes are indexed by a hashing of the substroke
        /// guids contained in the group.
        /// </summary>
        /// <param name="combinations">different lists of substrokes to compute bboxes for</param>
        /// <returns>Dictionary of group hash to group bounding box</returns>
        private Dictionary<int, System.Drawing.Rectangle> GetGroupBoxes(Dictionary<int, List<Substroke>> combinations)
        {
            Dictionary<int, System.Drawing.Rectangle> boxes = new Dictionary<int, System.Drawing.Rectangle>();
            foreach (KeyValuePair<int, List<Substroke>> kv in combinations)
            {
                List<Substroke> strokes = kv.Value;
                System.Drawing.Rectangle box = Utilities.Compute.BoundingBox(strokes.ToArray());
                boxes.Add(kv.Key, box);
            }

            return boxes;
        }

        /// <summary>
        /// Gets a hash code for a list of strokes by XORing (^)
        /// their GUIDs
        /// </summary>
        /// <param name="strokes">Strokes in group</param>
        /// <returns>Hash code for this group</returns>
        private int GetHashStrokes(List<Substroke> strokes)
        {
            int hash = 1;
            foreach (Substroke s in strokes)
                hash = hash ^ s.Id.GetHashCode();

            return hash;
        }

        /// <summary>
        /// Gets a hash code for bounding box, unique in position and size
        /// </summary>
        /// <param name="box"></param>
        /// <returns></returns>
        private int GetHashBbox(System.Drawing.Rectangle box)
        {
            return box.X ^ box.Y ^ box.Width ^ box.Height;
        }

        /// <summary>
        /// Populates a template's SymbolInfo using a matching template
        /// </summary>
        /// <param name="bestMatch">Best result from recognition</param>
        private void FillInSymbolInfo(ImageTemplateResult bestMatch)
        {
            // Fill out completeness of symbol
            if (bestMatch.Errors.Count == 0)
                m_SymbolInfo.Completeness = SymbolCompleteness.Complete;
            else
            {
                bool extra = false;
                bool missing = false;
                foreach (ImageMatchError error in bestMatch.Errors)
                {
                    if (error.Type == ErrorType.Missing)
                        missing = true;
                    else if (error.Type == ErrorType.Extra)
                        extra = true;
                }
                if (missing && extra)
                    m_SymbolInfo.Completeness = SymbolCompleteness.Combo;
                else if (missing)
                    m_SymbolInfo.Completeness = SymbolCompleteness.Partial;
                else if (extra)
                    m_SymbolInfo.Completeness = SymbolCompleteness.HasExtra;
                else
                    m_SymbolInfo.Completeness = SymbolCompleteness.Complete;
            }

            // Get the symbol's class and type
            m_SymbolInfo.SymbolClass = General.GetClass(bestMatch.Name);
            m_SymbolInfo.SymbolType = bestMatch.Name;
        }

        /// <summary>
        /// Creates a map of spatial relationships between all strokes in the shape
        /// </summary>
        /// <param name="completeShape"></param>
        /// <returns></returns>
        private Dictionary<Substroke, Dictionary<string, List<Substroke>>> GetLocationRelations(Shape completeShape)
        {
            List<Substroke> strokes = completeShape.SubstrokesL;

            Dictionary<Substroke, Dictionary<string, List<Substroke>>> map = 
                new Dictionary<Substroke, Dictionary<string, List<Substroke>>>(strokes.Count);

            Dictionary<Substroke, System.Drawing.Rectangle> strokeBoxes = new Dictionary<Substroke, System.Drawing.Rectangle>(strokes.Count);

            foreach (Substroke s in strokes)
            {
                System.Drawing.Rectangle box = Utilities.Compute.BoundingBox(s.PointsAsSysPoints);
                strokeBoxes.Add(s, box);

                Dictionary<string, List<Substroke>> direction = new Dictionary<string, List<Substroke>>(4);
                direction.Add("Left", new List<Substroke>());
                direction.Add("Right", new List<Substroke>());
                direction.Add("Above", new List<Substroke>());
                direction.Add("Below", new List<Substroke>());
                map.Add(s, direction);
            }

            for (int i = 0; i < strokes.Count; i++)
            {
                for (int j = i + 1; j < strokes.Count; j++)
                {
                    if (IsLeftOf(strokeBoxes[strokes[i]], strokeBoxes[strokes[j]]))
                    {
                        map[strokes[i]]["Right"].Add(strokes[j]);
                        map[strokes[j]]["Left"].Add(strokes[i]);
                    }
                    else if (IsRightOf(strokeBoxes[strokes[i]], strokeBoxes[strokes[j]]))
                    {
                        map[strokes[i]]["Left"].Add(strokes[j]);
                        map[strokes[j]]["Right"].Add(strokes[i]);
                    }

                    if (IsAbove(strokeBoxes[strokes[i]], strokeBoxes[strokes[j]]))
                    {
                        map[strokes[i]]["Above"].Add(strokes[j]);
                        map[strokes[j]]["Below"].Add(strokes[i]);
                    }
                    else if (IsBelow(strokeBoxes[strokes[i]], strokeBoxes[strokes[j]]))
                    {
                        map[strokes[i]]["Below"].Add(strokes[j]);
                        map[strokes[j]]["Above"].Add(strokes[i]);
                    }
                }
            }

            return map;
        }

        /// <summary>
        /// Is box2 below box1?
        /// </summary>
        /// <param name="box1"></param>
        /// <param name="box2"></param>
        /// <returns></returns>
        private bool IsBelow(System.Drawing.Rectangle box1, System.Drawing.Rectangle box2)
        {
            if (box1.Height > box2.Height)
            {
                if (box2.Top > box1.Bottom - box1.Height * MAX_BBOX_INCLUSION)
                    return true;
                else
                    return false;
            }
            else
            {
                if (box1.Bottom < box2.Top + box2.Height * MAX_BBOX_INCLUSION)
                    return true;
                else
                    return false;
            }
        }

        /// <summary>
        /// Is box2 above box1?
        /// </summary>
        /// <param name="box1"></param>
        /// <param name="box2"></param>
        /// <returns></returns>
        private bool IsAbove(System.Drawing.Rectangle box1, System.Drawing.Rectangle box2)
        {
            if (box2.Height > box1.Height)
            {
                if (box1.Top > box2.Bottom - box2.Height * MAX_BBOX_INCLUSION)
                    return true;
                else
                    return false;
            }
            else
            {
                if (box2.Bottom < box1.Top + box1.Height * MAX_BBOX_INCLUSION)
                    return true;
                else
                    return false;
            }
        }

        /// <summary>
        /// Is Box2 right of Box1
        /// </summary>
        /// <param name="box1"></param>
        /// <param name="box2"></param>
        /// <returns></returns>
        private bool IsRightOf(System.Drawing.Rectangle box1, System.Drawing.Rectangle box2)
        {
            if (box2.Width > box1.Width)
            {
                if (box1.Right < box2.Left + box2.Width * MAX_BBOX_INCLUSION)
                    return true;
                else
                    return false;
            }
            else
            {
                if (box2.Left > box1.Right - box1.Width * MAX_BBOX_INCLUSION)
                    return true;
                else
                    return false;
            }
        }

        /// <summary>
        /// Is Box2 left of Box1?
        /// </summary>
        /// <param name="box1"></param>
        /// <param name="box2"></param>
        /// <returns></returns>
        private bool IsLeftOf(System.Drawing.Rectangle box1, System.Drawing.Rectangle box2)
        {
            if (box2.Width > box1.Width)
            {
                if (box1.Right < box2.Left + box2.Width * MAX_BBOX_INCLUSION)
                    return true;
                else
                    return false;
            }
            else
            {
                if (box2.Left > box1.Right - box1.Width * MAX_BBOX_INCLUSION)
                    return true;
                else
                    return false;
            }
        }

        #endregion


        #region Getters

        /// <summary>
        /// Name of the template (e.g. AND, OR, Arrow)
        /// </summary>
        public string Name
        {
            get { return m_SymbolInfo.SymbolType; }
        }

        /// <summary>
        /// Shape parts found in this template (e.g. BackArc, Bubble, TopLine)
        /// </summary>
        public List<string> Parts
        {
            get { return new List<string>(m_ShapeParts.Values); }
        }

        /// <summary>
        /// Substrokes in the template
        /// </summary>
        public List<Substroke> Strokes
        {
            get { return new List<Substroke>(m_ShapeParts.Keys); }
        }

        /// <summary>
        /// Get's the Single-Stroke Recognition result (or known label)
        /// of the given stroke, returning "Unknown" if the substroke is
        /// not found in the template.
        /// </summary>
        /// <param name="stroke"></param>
        /// <returns></returns>
        public string GetShapePart(Substroke stroke)
        {
            if (m_ShapeParts.ContainsKey(stroke))
                return m_ShapeParts[stroke];
            else
                return "Unknown";
        }

        /// <summary>
        /// Indicates whether the template is known (i.e. whether it was
        /// generated by hand-labeled data).
        /// </summary>
        public bool IsFullyKnown
        {
            get { return m_IsFullyKnown; }
            set { m_IsFullyKnown = true; }
        }

        #endregion
    }
}
