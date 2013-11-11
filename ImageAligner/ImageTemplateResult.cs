/*
 * File: ImageTemplateResult.cs
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

namespace ImageAligner
{
    /// <summary>
    /// Type of error in the template matching
    /// </summary>
    public enum TemplateMatchingError { ExtraStroke, MissingStroke };

    /// <summary>
    /// Confidence levels in matching the templates
    /// </summary>
    public enum TemplateMatchingConfidence { High, Medium, Low, None };

    /// <summary>
    /// This object contains the results from matching two templates to each other. 
    /// The score indicates the likelihood that the match is correct, higher is better.
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
    [DebuggerDisplay("{m_Template.Name}: {m_Score} (Confidence: {m_Confidence})")]
    [Serializable]
    public class ImageTemplateResult
    {
        #region Constants / Scaling Parameters

        /// <summary>
        /// Scaling factor to keep scores in a nice range
        /// </summary>
        const double HAUS_SCALE = 1.0 / 48.0;

        /// <summary>
        /// Scaling factor to keep scores in a nice range
        /// </summary>
        const double MODHAUS_SCALE = 2.0 / 48.0;

        /// <summary>
        /// Penalty assessed to the score for each unmapped stroke
        /// </summary>
        const double PENALTY = 0.15;

        #endregion


        #region Member Variables

        /// <summary>
        /// Reference to the 'unknown' template
        /// </summary>
        ImageTemplate m_Unknown;

        /// <summary>
        /// Reference to the 'known' template that this result is for
        /// </summary>
        ImageTemplate m_Template;

        /// <summary>
        /// Mapping between strokes
        /// </summary>
        StrokeMapping m_Map;

        /// <summary>
        /// Individual results for each mapping pair
        /// </summary>
        Dictionary<ImageSymbol, ImageSymbolResult> m_SymbolResults;

        /// <summary>
        /// Combined Hausdorff distance for all symbols.
        /// Take the max of all individual Hausdorff distances.
        /// </summary>
        double m_HausdorffDistance;

        /// <summary>
        /// Combined Modified Hausdorff distance for all symbols.
        /// Take the average of all individual Modified Hausdorff distances.
        /// </summary>
        double m_ModifiedHausdorffDistance;

        /// <summary>
        /// Overall score (includes penalties for unmapped)
        /// </summary>
        double m_Score = -1.0;

        /// <summary>
        /// Image score (does NOT include penalties for unmapped)
        /// </summary>
        double m_ImageOnlyScore = -1.0;

        /// <summary>
        /// List of all errors between mappings
        /// </summary>
        List<ImageMatchError> m_Errors;

        /// <summary>
        /// Confidence in this mapping (based solely on overall score)
        /// </summary>
        TemplateMatchingConfidence m_Confidence;

        #endregion


        #region Constructors

        /// <summary>
        /// Constructor - Calculates scores after initializing
        /// </summary>
        /// <param name="unknown">ImageTemplate for the unknown (recognized) shape</param>
        /// <param name="template">ImageTemplate for the corresponding 'known' template</param>
        /// <param name="map">Stroke mapping between templates</param>
        /// <param name="results">Individual results from the stroke-mapping pairs</param>
        public ImageTemplateResult(ImageTemplate unknown, ImageTemplate template,
            StrokeMapping map, Dictionary<ImageSymbol, ImageSymbolResult> results)
        {
            m_Unknown = unknown;
            m_Template = template;
            m_Map = map;
            m_SymbolResults = results;
            m_Errors = new List<ImageMatchError>();

            CalculateScores();
        }

        #endregion


        #region Calculation Functions

        /// <summary>
        /// Takes the individual results from each stroke-pairing and 
        /// combines them to form composite scores.
        /// </summary>
        private void CalculateScores()
        {
            // Combine individual Haussdorf and Modified Haussdorf distances
            double hauss = 0.0;
            double modHauss = 0.0;
            foreach (ImageSymbolResult result in m_SymbolResults.Values)
            {
                hauss = Math.Max(hauss, result.Hausdorff);
                modHauss += result.ModifiedHausdorff;
            }

            if (m_SymbolResults.Count > 0)
                modHauss /= m_SymbolResults.Count;
            else
            {
                hauss = double.PositiveInfinity;
                modHauss = double.PositiveInfinity;
            }
            m_HausdorffDistance = hauss;
            m_ModifiedHausdorffDistance = modHauss;

            // Identify Errors
            int penaltyErrorCount = 0;
            foreach (Substroke s in m_Unknown.Strokes)
                if (!m_Map.ContainsKey(s))
                {
                    ErrorDetail detail = GetDetail(s, m_Unknown);
                    ErrorSeverity severity = GetSeverity(s, m_Unknown);
                    ImageMatchError error = new ImageMatchError(ErrorType.Extra, detail, severity, s);
                    m_Errors.Add(error);
                    penaltyErrorCount++;
                }

            foreach (Substroke s in m_Template.Strokes)
                if (!m_Map.ContainsValue(s))
                {
                    ErrorDetail detail = GetDetail(s, m_Template);
                    ErrorSeverity severity = GetSeverity(s, m_Template);
                    ImageMatchError error = new ImageMatchError(ErrorType.Missing, detail, severity);
                    m_Errors.Add(error);
                }

            double penalty = penaltyErrorCount * PENALTY;

            // Compute Scores
            m_ImageOnlyScore = 1.0 / (1.0 + m_HausdorffDistance * HAUS_SCALE + m_ModifiedHausdorffDistance * MODHAUS_SCALE); 
            m_Score = m_ImageOnlyScore - penalty;
            
            // Determine Confidence - contrived values
            if (m_Score > 0.83)
                m_Confidence = TemplateMatchingConfidence.High;
            else if (m_Score > 0.7)
                m_Confidence = TemplateMatchingConfidence.Medium;
            else if (m_Score > 0.57)
                m_Confidence = TemplateMatchingConfidence.Low;
            else
                m_Confidence = TemplateMatchingConfidence.None;
        }

        /// <summary>
        /// Determine the severity of an error based on size and location
        /// </summary>
        /// <param name="stroke">Stroke with associated error to determine severity for</param>
        /// <param name="template">ImateTemplate - either 'known' or 'unknown'</param>
        /// <returns>Severity of error</returns>
        private ErrorSeverity GetSeverity(Substroke stroke, ImageTemplate template)
        {
            // If there problem is just that the unknown is missing a 
            // junk or touchup stroke, the severity is low
            if (template.IsFullyKnown)
            {
                string part = template.GetShapePart(stroke);
                if (part == "Junk" || part == "TouchUp")
                    return ErrorSeverity.Low;
            }

            // Determine the bounding box for the stroke in question, as well
            // as the shape when excluding the stroke in question
            System.Drawing.Rectangle bbox = Utilities.Compute.BoundingBox(stroke.PointsAsSysPoints);
            List<Substroke> strokes = new List<Substroke>(template.Strokes);
            strokes.Remove(stroke);
            System.Drawing.Rectangle shapebox = Utilities.Compute.BoundingBox(strokes.ToArray());

            // Get the proportional size and the amount of overlap between the boxes
            double proportion = Utilities.Compute.AreaComparison_Percent_R1toR2(bbox, shapebox);
            double percentInside = Utilities.Compute.PercentRectangle1InsideRectangle2(bbox, shapebox);

            // Assign severities
            if (percentInside > 0.9)
            {
                if (proportion < 0.5)
                    return ErrorSeverity.Low;
                else
                    return ErrorSeverity.Medium;
            }
            else if (percentInside > 0.1)
            {
                if (proportion < 0.3)
                    return ErrorSeverity.Medium;
                else
                    return ErrorSeverity.High;
            }
            else
            {
                if (proportion > 0.1)
                    return ErrorSeverity.High;
                else
                    return ErrorSeverity.Medium;
            }
        }

        /// <summary>
        /// Gets the name of the shape part that is missing or extra
        /// </summary>
        /// <param name="stroke">Stroke to get part name of</param>
        /// <param name="template">ImageTemplate containing the stroke</param>
        /// <returns>Shape part for the stroke</returns>
        private ErrorDetail GetDetail(Substroke stroke, ImageTemplate template)
        {
            string type = template.GetShapePart(stroke);
            return ImageTemplateResult.GetDetail(type);
        }

        public static ErrorDetail GetDetail(string type)
        {
            switch (type)
            {
                case "Unknown":
                    return ErrorDetail.Unknown;
                case "BackArc":
                    return ErrorDetail.BackArc;
                case "BackLine":
                    return ErrorDetail.BackLine;
                case "BottomArc":
                    return ErrorDetail.BottomArc;
                case "BottomLine":
                    return ErrorDetail.BottomLine;
                case "Bubble":
                    return ErrorDetail.Bubble;
                case "FrontArc":
                    return ErrorDetail.FrontArc;
                case "GreaterThan":
                    return ErrorDetail.GreaterThan;
                case "Junk":
                    return ErrorDetail.Junk;
                case "Label":
                    return ErrorDetail.Label;
                case "LabelBoxOther":
                    return ErrorDetail.LabelBoxOther;
                case "Not_Hat":
                    return ErrorDetail.Not_Hat;
                case "Not_V":
                    return ErrorDetail.Not_V;
                case "TopArc":
                    return ErrorDetail.TopArc;
                case "TopLine":
                    return ErrorDetail.TopLine;
                case "TouchUp":
                    return ErrorDetail.TouchUp;
                case "Triangle":
                    return ErrorDetail.Triangle;
                case "Wire":
                    return ErrorDetail.Wire;
                default:
                    return ErrorDetail.Unknown;
            }
        }

        #endregion


        #region Getters

        /// <summary>
        /// Overall score - accounts for mapping errors
        /// </summary>
        public double Score
        {
            get { return m_Score; }
        }

        /// <summary>
        /// Score for the image comparison only - ignores mapping errors
        /// </summary>
        public double ImageOnlyScore
        {
            get { return m_ImageOnlyScore; }
        }

        /// <summary>
        /// Gets the template's name (type) (e.g. AND, OR)
        /// </summary>
        public string Name
        {
            get 
            {
                if (m_Template != null)
                    return m_Template.Name;
                else
                    return "Unknown";
            }
        }

        /// <summary>
        /// Gets the errors found in this matching
        /// </summary>
        public List<ImageMatchError> Errors
        {
            get { return m_Errors; }
        }

        /// <summary>
        /// Provides a way to remove an error from the list, 
        /// indicating that it has been resolve.
        /// </summary>
        /// <param name="error"></param>
        public void RemoveError(ImageMatchError error)
        {
            if (m_Errors.Contains(error))
            {
                m_Errors.Remove(error);
                m_Score += PENALTY;
            }
        }

        /// <summary>
        /// Gets the confidence in the matching's correctness
        /// </summary>
        public TemplateMatchingConfidence Confidence
        {
            get { return m_Confidence; }
        }

        /// <summary>
        /// Gets the mapping of strokes between the two templates
        /// </summary>
        public StrokeMapping Map
        {
            get { return m_Map; }
        }

        public string GetShapePart(Substroke stroke)
        {
            return m_Template.GetShapePart(stroke);
        }

        public List<Substroke> SubstrokesUsedInMatch
        {
            get
            {
                List<Substroke> strokes = new List<Substroke>();

                foreach (Substroke stroke in m_Map.StrokeValues)
                    strokes.Add(stroke);

                return strokes;
            }
        }

        public List<Substroke> SubstrokesNOTUsedInMatch
        {
            get
            {
                List<Substroke> strokes = new List<Substroke>();

                foreach (Substroke stroke in m_Map.UnMappedStrokeValues)
                    strokes.Add(stroke);

                return strokes;
            }
        }

        #endregion
    }
}
