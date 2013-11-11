using System;
using System.Collections.Generic;
using System.Runtime.Serialization;
using System.Text;
using System.Diagnostics;
using Sketch;

namespace ImageAligner
{
    public enum GroupingCorrectness { CompleteParentShape, CompleteSubShape,
        IncorrectParentShapeMissing, IncorrectParentShapeExtra, IncorrectParentShapeMissingAndExtra,
        IncorrectSubShapeMissing, IncorrectSubShapeExtra, IncorrectSubShapeMissingAndExtra,
        NoShapeOfCorrectType };
    public enum ResultCompleteness { Correct, ExtraErrors, MissingErrors, BothMissingAndExtraErrors, NoResult };
    public enum ResultConfidence { CorrectHighMidConfidence, CorrectLowConfidence,
        IncorrectHighMidConfidence, IncorrectLowConfidence, NoResult };
    public enum ResultType { Correct, InTop3, Incorrect, NoResult };
    enum completion { Complete, Partial, None };


    [Serializable]
    public class TestSketchAccuracy
    {
        Sketch.Sketch m_Sketch;
        Dictionary<Shape, List<ImageTemplateResult>> m_Shapes;
        List<ShapeResult> m_Results;

        Dictionary<Shape, Shape> m_LinkShapes;

        public TestSketchAccuracy(Sketch.Sketch sketch, Dictionary<Shape, List<ImageTemplateResult>> shapes)
        {
            m_LinkShapes = new Dictionary<Shape, Shape>();
            sketch = Utilities.General.ReOrderParentShapes(sketch);
            m_Sketch = Utilities.General.LinkShapes(sketch);
            m_Shapes = shapes;
            m_Results = new List<ShapeResult>();

            CalculateAccuracy2();
        }

        private void CalculateAccuracy2()
        {
            foreach (Shape shape in m_Shapes.Keys)
            {
                Shape parentOfBest;
                Shape bestMatchingShape = FindBestMatchingShape(shape, m_Sketch, out parentOfBest);

                ShapeResult result = new ShapeResult(shape, bestMatchingShape, parentOfBest, m_Shapes[shape]);
                m_Results.Add(result);
            }
        }

        private Shape FindBestMatchingShape(Shape shape, Sketch.Sketch sketch, out Shape ParentOfBest)
        {
            List<Shape> parents = new List<Shape>();
            Dictionary<Shape, completion> completeness = new Dictionary<Shape, completion>();
            Dictionary<Shape, int> numStrokesInShape = new Dictionary<Shape, int>();

            foreach (Substroke stroke in shape.SubstrokesL)
            {
                Substroke realStroke = sketch.GetSubstroke(stroke.Id);
                foreach (Shape parent in realStroke.ParentShapes)
                {
                    if (!completeness.ContainsKey(parent) && parent.ShapesL.Count > 0 && Utilities.General.IsGate(parent))
                    {
                        completion comp = FindCompleteness(parent, shape);
                        completeness.Add(parent, comp);
                        int count = NumStrokesInParent(shape, parent);
                        numStrokesInShape.Add(parent, count);
                        parents.Add(parent);
                    }
                }
            }

            if (parents.Count == 0)
            {
                foreach (Substroke stroke in shape.SubstrokesL)
                {
                    Substroke realStroke = sketch.GetSubstroke(stroke.Id);
                    foreach (Shape parent in realStroke.ParentShapes)
                    {
                        if (!completeness.ContainsKey(parent))
                        {
                            completion comp = FindCompleteness(parent, shape);
                            completeness.Add(parent, comp);
                            int count = NumStrokesInParent(shape, parent);
                            numStrokesInShape.Add(parent, count);
                            parents.Add(parent);
                        }
                    }
                }
            }

            Shape best = null;
            int bestCount = 0;
            completion bestComplete = completion.None;
            bool bestIsSubshape = false;
            foreach (Shape parent in parents)
            {
                int num = numStrokesInShape[parent];
                completion complete = completeness[parent];
                bool isParentShape = IsShapeNull(parent.ParentShape);

                if (num > bestCount)
                {
                    bestCount = num;
                    bestComplete = complete;
                    bestIsSubshape = !isParentShape;
                    best = parent;
                }
                else if (num == bestCount)
                {
                    if (bestComplete == completion.Complete)
                    {
                        continue;
                    }
                    else if (bestComplete == completion.Partial)
                    {
                        if (complete == completion.Complete)
                        {
                            bestCount = num;
                            bestComplete = complete;
                            bestIsSubshape = !isParentShape;
                            best = parent;
                        }
                        else if (complete == completion.Partial)
                        {
                            if (!bestIsSubshape && !isParentShape)
                            {
                                bestCount = num;
                                bestComplete = complete;
                                bestIsSubshape = !isParentShape;
                                best = parent;
                            }
                        }
                    }
                    else if (bestComplete == completion.None)
                    {
                        bestCount = num;
                        bestComplete = complete;
                        bestIsSubshape = !isParentShape;
                        best = parent;
                    }
                }
            }

            ParentOfBest = null;
            if (best != null && !IsShapeNull(best.ParentShape))
                ParentOfBest = best.ParentShape;

            return best;
        }

        private bool IsShapeNull(Shape shape)
        {
            if (shape.SubstrokesL.Count == 0)
                return true;
            else
                return false;
        }

        private completion FindCompleteness(Shape parent, Shape shape)
        {
            bool foundAny = false;
            bool missingAny = false;

            foreach (Substroke stroke in parent.SubstrokesL)
            {
                bool found = Utilities.General.ShapeContainsSubstroke(shape, stroke);
                if (found)
                    foundAny = true;
                else
                    missingAny = true;
            }

            if (!foundAny)
                return completion.None;
            else if (missingAny)
                return completion.Partial;
            else
                return completion.Complete;
        }

        public List<ShapeResult> Results
        {
            get { return m_Results; }
        }


        #region old

        
        /*private void CalculateAccuracy()
        {
            List<Shape> SketchParentShapes = new List<Shape>();
            List<Shape> SketchSubShapes = new List<Shape>();
            List<Shape> SketchLeafShapes = new List<Shape>();
            foreach (Shape shape in m_Sketch.ShapesL)
            {
                if (!Utilities.General.IsGate(shape) || Utilities.General.leafLabels.Contains(shape.Type))
                    continue;

                if (shape == shape.SubstrokesL[0].ParentShapes[0])
                    SketchParentShapes.Add(shape);
                else if (shape == shape.SubstrokesL[0].ParentShapes[shape.SubstrokesL[0].ParentShapes.Count - 1])
                    SketchLeafShapes.Add(shape);
                else
                    SketchSubShapes.Add(shape);
            }

            foreach (Shape shape in m_Shapes.Keys)
            {
                ShapeCorrectness correctness = new ShapeCorrectness();

                #region Find Grouping Correctness

                Shape PerfectParent = null;
                Shape PerfectSubShape = null;
                Shape BestMatchingParentShape = null;
                Shape BestMatchingSubShape = null;
                List<Substroke> missingStrokes = new List<Substroke>();
                List<Substroke> extraStrokes = new List<Substroke>();

                Shape matchingShape;
                bool found = FindEquivalentShape(shape, m_Sketch, out matchingShape);

                if (found)
                {
                    #region Perfect matching shape
                    // We have a matching shape from the labeled data
                    if (matchingShape == matchingShape.Substrokes[0].ParentShapes[0])
                    {
                        // The matching shape is the highest parent (e.g. matching=NOR, NOR contains OR)
                        correctness.Grouping_Correctness = GroupingCorrectness.CompleteParentShape;
                        PerfectParent = matchingShape;
                    }
                    else if (matchingShape == matchingShape.Substrokes[0].ParentShapes[matchingShape.Substrokes[0].ParentShapes.Count - 1])
                    {
                        // Should never get here
                        // The matching shape is the leaf shape(e.g. BackArc)
                        //correctness.Grouping_Correctness = GroupingCorrectness.SingleStrokeGroup;
                        //PerfectLeafShape = matchingShape;
                    }
                    else
                    {
                        // The matching shape is a middle shape (e.g. matching=OR, this OR is contained in a NOR)
                        correctness.Grouping_Correctness = GroupingCorrectness.CompleteSubShape;
                        PerfectSubShape = matchingShape;
                    }
                    #endregion
                }
                else
                {
                    // Need to add check for case such as NOTBUBBLE (parent) has
                    // been grouped with strokes from another (incomplete parent) gate
                    #region Check for Parent Shape with extra strokes
                    List<Shape> completelyContainedParents = new List<Shape>();
                    foreach (Shape parent in SketchParentShapes)
                        if (Shape1CompletelyContainedInShape2(parent, shape))
                            completelyContainedParents.Add(parent);

                    bool done = false;
                    if (completelyContainedParents.Count == 1)
                    {
                        BestMatchingParentShape = completelyContainedParents[0];
                        correctness.Grouping_Correctness = GroupingCorrectness.IncorrectParentShapeExtra;

                        List<Guid> ids = new List<Guid>(BestMatchingParentShape.SubstrokesL.Count);
                        foreach (Substroke stroke in BestMatchingParentShape.SubstrokesL)
                            ids.Add(stroke.Id);

                        foreach (Substroke stroke in shape.SubstrokesL)
                            if (!ids.Contains(stroke.Id))
                                extraStrokes.Add(stroke);

                        done = true;
                    }
                    else if (completelyContainedParents.Count > 1)
                    {
                        int highest = 0;
                        foreach (Shape parent in completelyContainedParents)
                        {
                            if (parent.SubstrokesL.Count > highest)
                            {
                                highest = parent.SubstrokesL.Count;
                                BestMatchingParentShape = parent;
                                correctness.Grouping_Correctness = GroupingCorrectness.IncorrectParentShapeExtra;
                            }
                        }

                        List<Guid> ids = new List<Guid>(BestMatchingParentShape.SubstrokesL.Count);
                        foreach (Substroke stroke in BestMatchingParentShape.SubstrokesL)
                            ids.Add(stroke.Id);

                        foreach (Substroke stroke in shape.SubstrokesL)
                            if (!ids.Contains(stroke.Id))
                                extraStrokes.Add(stroke);

                        done = true;
                    }
                    #endregion

                    #region Check for missing strokes
                    if (!done)
                    {
                        int highest = 0;
                        foreach (Shape parent in SketchParentShapes)
                        {
                            int numFound = NumStrokesInParent(shape, parent);
                            if (numFound > highest)
                            {
                                highest = numFound;
                                BestMatchingParentShape = parent;
                            }
                        }

                        if (highest == 0 || BestMatchingParentShape == null)
                        {
                            // Found a group that only contains wires or labels
                            correctness.Grouping_Correctness = GroupingCorrectness.NoShapeOfCorrectType;
                        }
                        else
                        {
                            // Check whether this is a complete subshape with extra strokes not part
                            // of the parent
                            bool subFound = false;
                            int numStrokesInSub = 0;
                            foreach (Shape subShape in SketchSubShapes)
                            {
                                if (Shape1CompletelyContainedInShape2(subShape, shape))
                                {
                                    List<Substroke> otherStrokes = new List<Substroke>();
                                    List<Guid> ids = new List<Guid>(shape.SubstrokesL.Count);
                                    foreach (Substroke stroke in subShape.SubstrokesL)
                                        ids.Add(stroke.Id);

                                    foreach (Substroke stroke in shape.SubstrokesL)
                                        if (!ids.Contains(stroke.Id))
                                            otherStrokes.Add(stroke);

                                    Shape parentOfSubShape = subShape.SubstrokesL[0].ParentShapes[0];
                                    List<Guid> idsParent = new List<Guid>(parentOfSubShape.SubstrokesL.Count);
                                    foreach (Substroke stroke in parentOfSubShape.SubstrokesL)
                                        idsParent.Add(stroke.Id);

                                    bool foundGoodSub = true;
                                    foreach (Substroke stroke in otherStrokes)
                                    {
                                        if (idsParent.Contains(stroke.Id))
                                            foundGoodSub = false;
                                    }

                                    if (foundGoodSub && subShape.SubstrokesL.Count > numStrokesInSub)
                                    {
                                        numStrokesInSub = subShape.SubstrokesL.Count;
                                        BestMatchingParentShape = null;
                                        subFound = true;
                                        BestMatchingSubShape = subShape;
                                    }

                                }
                            }

                            if (!subFound)
                            {
                                List<Guid> ids = new List<Guid>(BestMatchingParentShape.SubstrokesL.Count);
                                foreach (Substroke stroke in BestMatchingParentShape.SubstrokesL)
                                    ids.Add(stroke.Id);

                                foreach (Substroke stroke in shape.SubstrokesL)
                                    if (!ids.Contains(stroke.Id))
                                        extraStrokes.Add(stroke);

                                List<Guid> ids2 = new List<Guid>(shape.SubstrokesL.Count);
                                foreach (Substroke stroke in shape.SubstrokesL)
                                    ids2.Add(stroke.Id);

                                foreach (Substroke stroke in BestMatchingParentShape.SubstrokesL)
                                    if (!ids2.Contains(stroke.Id))
                                        missingStrokes.Add(stroke);

                                if (missingStrokes.Count > 0 && extraStrokes.Count > 0)
                                    correctness.Grouping_Correctness = GroupingCorrectness.IncorrectParentShapeMissingAndExtra;
                                else if (missingStrokes.Count > 0)
                                    correctness.Grouping_Correctness = GroupingCorrectness.IncorrectParentShapeMissing;
                                else if (extraStrokes.Count > 0)
                                    correctness.Grouping_Correctness = GroupingCorrectness.IncorrectParentShapeExtra;
                            }
                            else
                            {
                                List<Guid> ids = new List<Guid>(BestMatchingSubShape.SubstrokesL.Count);
                                foreach (Substroke stroke in BestMatchingSubShape.SubstrokesL)
                                    ids.Add(stroke.Id);

                                foreach (Substroke stroke in shape.SubstrokesL)
                                    if (!ids.Contains(stroke.Id))
                                        extraStrokes.Add(stroke);

                                List<Guid> ids2 = new List<Guid>(shape.SubstrokesL.Count);
                                foreach (Substroke stroke in shape.SubstrokesL)
                                    ids2.Add(stroke.Id);

                                foreach (Substroke stroke in BestMatchingSubShape.SubstrokesL)
                                    if (!ids2.Contains(stroke.Id))
                                        missingStrokes.Add(stroke);

                                if (missingStrokes.Count > 0 && extraStrokes.Count > 0)
                                    correctness.Grouping_Correctness = GroupingCorrectness.IncorrectSubShapeMissingAndExtra;
                                else if (missingStrokes.Count > 0)
                                    correctness.Grouping_Correctness = GroupingCorrectness.IncorrectSubShapeMissing;
                                else if (extraStrokes.Count > 0)
                                    correctness.Grouping_Correctness = GroupingCorrectness.IncorrectSubShapeExtra;
                            }
                        }
                    }*/

                    /*
                    // Get all the parent shapes
                    List<Shape> parentShapes = new List<Shape>();
                    List<Shape> subShapes = new List<Shape>();
                    List<Shape> leafShapes = new List<Shape>();
                    foreach (Substroke stroke in shape.Substrokes)
                    {
                        Substroke labeledStroke = m_Sketch.GetSubstroke(stroke.Id);
                        for (int i = 0; i < labeledStroke.ParentShapes.Count; i++)
                        {
                            Shape current = labeledStroke.ParentShapes[i];
                            if (i == 0)
                                parentShapes.Add(current);
                            else if (i == labeledStroke.ParentShapes.Count - 1)
                                leafShapes.Add(current);
                            else
                                subShapes.Add(current);
                        }
                    }

                    Dictionary<Shape, int> parentCount = new Dictionary<Shape, int>();
                    // Find best matching gate shape
                    foreach (Shape parent in parentShapes)
                    {
                        if (!Utilities.General.IsGate(parent.Type)) continue;

                        if (!parentCount.ContainsKey(parent))
                            parentCount.Add(parent, 0);

                        parentCount[parent]++;
                    }
                    Shape bestParent = null;
                    

                    if (parentCount.Count == 0)
                    {
                        // No matching gate shape...
                    }
                    else if (parentCount.Count > 1)
                    {
                        // Multiple shapes
                        int bestScore = 0;
                        foreach (KeyValuePair<Shape, int> kvp in parentCount)
                        {
                            if (kvp.Value > bestScore)
                            {
                                bestScore = kvp.Value;
                                bestParent = kvp.Key;
                            }
                            else if (kvp.Value == bestScore)
                            {
                                // What should I do if there are equal numbers?
                            }
                        }
                    }
                    else
                    {
                        // Single parent shape
                        bestParent = new List<Shape>(parentCount.Keys)[0];
                    }

                    bestMatchingParentShape = bestParent;

                    // Now that we've found the best matching parent...
                    if (bestParent != null)
                    {
                        // Any missing strokes in 'shape'?
                        foreach (Substroke stroke in bestParent.SubstrokesL)
                        {
                            if (!ShapeContainsSubstroke(shape, stroke))
                                missingStrokes.Add(stroke);
                        }

                        // Any extra strokes in 'shape'?
                        foreach (Substroke stroke in shape.SubstrokesL)
                        {
                            if (!ShapeContainsSubstroke(bestParent, stroke))
                                extraStrokes.Add(stroke);
                        }
                    }


                    if (missingStrokes.Count > 0)
                    {
                        List<Shape> parentSubs = new List<Shape>();
                        foreach (Substroke stroke in bestParent.SubstrokesL)
                        {
                            for (int i = 1; i < stroke.ParentShapes.Count; i++)
                            {
                                if (i != stroke.ParentShapes.Count - 1)
                                    parentSubs.Add(stroke.ParentShapes[i]);
                            }
                        }

                        int highest = 0;
                        Shape bestSub = null;
                        foreach (Shape pSub in parentSubs)
                        {
                            bool allFound = true;
                            foreach (Substroke stroke in pSub.SubstrokesL)
                            {
                                if (!ShapeContainsSubstroke(shape, stroke))
                                    allFound = false;
                            }

                            if (allFound && pSub.SubstrokesL.Count > highest)
                            {
                                highest = pSub.SubstrokesL.Count;
                                bestSub = pSub;
                            }
                        }


                        if (bestSub == null)
                        {
                            Dictionary<Shape, int> subShapeCount = new Dictionary<Shape, int>();
                            // Check if any subshapes are completely found in 'shape'
                            foreach (Shape subShape in subShapes)
                            {
                                if (!Utilities.General.IsGate(subShape.Type)) continue;

                                if (!subShapeCount.ContainsKey(subShape))
                                    subShapeCount.Add(subShape, 0);

                                subShapeCount[subShape]++;
                            }
                            Shape bestSubShape = null;


                        }
                    }
                    
                    
                }



                #region Result Completeness and Confidence

                List<ImageTemplateResult> results = m_Shapes[shape];
                if (results.Count == 0)
                {
                    correctness.Result_Completeness = ResultCompleteness.NoResult;
                    correctness.Result_Confidence = ResultConfidence.NoResult;
                    correctness.Result_Type = ResultType.NoResult;
                    continue;
                }

                ImageTemplateResult topResult = results[0];
                if (correctness.Grouping_Correctness == GroupingCorrectness.CompleteParentShape)
                {
                    if (topResult.Errors.Count == 0)
                        correctness.Result_Completeness = ResultCompleteness.Correct;
                    else
                        correctness.Result_Completeness = ResultCompleteness.ExtraErrors;

                    if (topResult.Confidence != TemplateMatchingConfidence.Low)
                        correctness.Result_Confidence = ResultConfidence.CorrectHighMidConfidence;
                    else
                        correctness.Result_Confidence = ResultConfidence.IncorrectLowConfidence;
                }
                else if (correctness.Grouping_Correctness == GroupingCorrectness.CompleteSubShape)
                {
                    if (topResult.Errors.Count == 0)
                        correctness.Result_Completeness = ResultCompleteness.Correct;
                    else
                        correctness.Result_Completeness = ResultCompleteness.ExtraErrors;

                    if (topResult.Confidence != TemplateMatchingConfidence.Low)
                        correctness.Result_Confidence = ResultConfidence.CorrectHighMidConfidence;
                    else
                        correctness.Result_Confidence = ResultConfidence.IncorrectLowConfidence;
                }
                else if (correctness.Grouping_Correctness == GroupingCorrectness.IncorrectParentShapeMissingAndExtra
                    || correctness.Grouping_Correctness == GroupingCorrectness.IncorrectParentShapeMissing
                    || correctness.Grouping_Correctness == GroupingCorrectness.IncorrectParentShapeExtra
                    || correctness.Grouping_Correctness == GroupingCorrectness.IncorrectSubShapeMissingAndExtra
                    || correctness.Grouping_Correctness == GroupingCorrectness.IncorrectSubShapeMissing
                    || correctness.Grouping_Correctness == GroupingCorrectness.IncorrectSubShapeExtra)
                {
                    if (topResult.Errors.Count == 0)
                        correctness.Result_Completeness = ResultCompleteness.MissingErrors;

                    List<ImageMatchError> shouldHaveErrors = new List<ImageMatchError>();
                    foreach (Substroke missing in missingStrokes)
                    {
                        ImageMatchError error = new ImageMatchError(ErrorType.Missing, ErrorDetail.BackArc, ErrorSeverity.High, missing);
                        shouldHaveErrors.Add(error);
                    }
                    foreach (Substroke extra in extraStrokes)
                    {
                        ImageMatchError error = new ImageMatchError(ErrorType.Extra, ErrorDetail.BackArc, ErrorSeverity.High, extra);
                        shouldHaveErrors.Add(error);
                    }

                    int extraErrors = 0;
                    int missingErrors = 0;

                    foreach (ImageMatchError error in topResult.Errors)
                    {
                        bool foundError = false;
                        foreach (ImageMatchError shouldError in shouldHaveErrors)
                        {
                            if (!foundError && error.Type == shouldError.Type && error.OffendingStroke == shouldError.OffendingStroke)
                                    foundError = true;
                        }
                        if (!foundError)
                            extraErrors++;
                    }

                    foreach (ImageMatchError shouldError in shouldHaveErrors)
                    {
                        bool foundError = false;
                        foreach (ImageMatchError error in topResult.Errors)
                        {
                            if (!foundError && error.Type == shouldError.Type && error.OffendingStroke == shouldError.OffendingStroke)
                                foundError = true;
                        }
                        if (!foundError)
                            missingErrors++;
                    }
                }
                else if (correctness.Grouping_Correctness == GroupingCorrectness.IncorrectSubShapeMissingAndExtra
                    || correctness.Grouping_Correctness == GroupingCorrectness.IncorrectSubShapeMissing
                    || correctness.Grouping_Correctness == GroupingCorrectness.IncorrectSubShapeExtra)
                {

                }
                else if (correctness.Grouping_Correctness == GroupingCorrectness.NoShapeOfCorrectType)
                {
                    correctness.Result_Completeness = ResultCompleteness.MissingErrors;
                    if (topResult.Confidence != TemplateMatchingConfidence.Low)
                        correctness.Result_Confidence = ResultConfidence.CorrectLowConfidence;
                    else
                        correctness.Result_Confidence = ResultConfidence.IncorrectHighMidConfidence;
                }

                #endregion

                #region Result Type

                string label = "";
                if (PerfectParent != null)
                    label = PerfectParent.Type;
                else if (PerfectSubShape != null)
                    label = PerfectSubShape.Type;
                else if (BestMatchingParentShape != null)
                    label = BestMatchingParentShape.Type;
                else if (BestMatchingSubShape != null)
                    label = BestMatchingSubShape.Type;
                
                if (topResult.Name == label)
                    correctness.Result_Type = ResultType.Correct;
                else
                {
                    for (int i = 1; i < results.Count && i < 3; i++)
                    {
                        if (results[i].Name == label)
                            correctness.Result_Type = ResultType.InTop3;
                    }
                    if (correctness.Result_Type != ResultType.InTop3)
                        correctness.Result_Type = ResultType.Incorrect;
                }

                #endregion

                m_Correctness.Add(shape, correctness);
            }
        }*/

                    #endregion

        #region helpers

        private int NumStrokesInParent(Shape shape, Shape parent)
        {
            int found = 0;
            List<Guid> ids = new List<Guid>(parent.SubstrokesL.Count);
            foreach (Substroke stroke in parent.SubstrokesL)
                ids.Add(stroke.Id);

            foreach (Substroke stroke in shape.SubstrokesL)
                if (ids.Contains(stroke.Id))
                    found++;

            return found;
        }

        private bool Shape1CompletelyContainedInShape2(Shape shape1, Shape shape2)
        {
            List<Guid> ids = new List<Guid>(shape2.SubstrokesL.Count);
            foreach (Substroke stroke in shape2.SubstrokesL)
                ids.Add(stroke.Id);

            foreach (Substroke stroke in shape1.SubstrokesL)
                if (!ids.Contains(stroke.Id))
                    return false;

            return true;
        }

        private bool FindEquivalentShape(Shape shape, Sketch.Sketch sketch, out Shape matchingShape)
        {
            foreach (Shape s in sketch.Shapes)
            {
                if (AreEquivalent(shape, s))
                {
                    if (!Utilities.General.leafLabels.Contains(s.Type))
                    {
                        matchingShape = s;
                        return true;
                    }
                }
            }

            matchingShape = null;
            return false;
        }

        private bool AreEquivalent(Shape shape, Shape s)
        {
            if (shape.Substrokes.Length != s.Substrokes.Length)
                return false;

            List<Guid> ids = new List<Guid>(s.SubstrokesL.Count);
            foreach (Substroke stroke in s.SubstrokesL)
                ids.Add(stroke.Id);

            foreach (Substroke stroke1 in shape.Substrokes)
            {
                if (!ids.Contains(stroke1.Id))
                    return false;
            }

            return true;
        }

        #endregion
    }

    [Serializable]
    [DebuggerDisplay("Expected: {m_BestMatchingShape.Label} w/ {m_ExpectedErrors.Count} Errors --> Actual: {TopResult.Name} w/ {TopResult.Errors.Count} Errors, Score = {TopResult.Score.ToString(\"#0.00\")} ({TopResult.Confidence.ToString()})")]
    public class ShapeResult : ISerializable
    {
        Shape m_GroupedShape;
        List<ImageTemplateResult> m_ImageResults;

        Shape m_BestMatchingShape;
        Shape m_BestMatchingShapeParent;
        List<ImageMatchError> m_ExpectedErrors;
        


        public ShapeResult(Shape grouped, Shape bestMatch, Shape bestMatchParent, List<ImageTemplateResult> results)
        {
            m_GroupedShape = grouped;
            m_BestMatchingShape = bestMatch;
            m_BestMatchingShapeParent = bestMatchParent;
            m_ImageResults = results;
            if (bestMatch != null)
                m_ExpectedErrors = FindExpectedErrors();         
        }

        public Shape GroupedShape
        {
            get { return m_GroupedShape; }
        }

        public Shape BestMatchingShape
        {
            get { return m_BestMatchingShape; }
        }

        public List<ImageMatchError> ExpectedErrors
        {
            get { return m_ExpectedErrors; }
        }

        public bool IsPerfect
        {
            get
            {
                if (this.ExpectedShapeIsParent && this.BothShapesComplete && this.IsCorrectType && this.TopResult.Confidence != TemplateMatchingConfidence.Low)
                    return true;
                else
                    return false;
            }
        }

        public bool IsHumanCorrect
        {
            get
            {
                if (this.IsCorrectType && this.DoesTopResultGivePerfectGuidance)
                    return true;
                else
                    return false;
            }
        }

        public bool BothShapesComplete
        {
            get
            {
                if (m_ExpectedErrors.Count != 0)
                    return false;
                if (TopResult == null)
                    return false;
                if (TopResult.Errors.Count != 0)
                    return false;

                return true;
            }
        }

        public ImageTemplateResult TopResult
        {
            get
            {
                if (m_ImageResults != null && m_ImageResults.Count > 0)
                    return m_ImageResults[0];
                else
                    return null;
            }
        }

        public ImageTemplateResult ResultNum2
        {
            get
            {
                if (m_ImageResults != null && m_ImageResults.Count > 1)
                    return m_ImageResults[1];
                else
                    return null;
            }
        }

        public ImageTemplateResult ResultNum3
        {
            get
            {
                if (m_ImageResults != null && m_ImageResults.Count > 2)
                    return m_ImageResults[2];
                else
                    return null;
            }
        }

        private List<ImageMatchError> FindExpectedErrors()
        {
            List<ImageMatchError> errors = new List<ImageMatchError>();
            foreach (Substroke stroke in m_GroupedShape.SubstrokesL)
            {
                if (!Utilities.General.ShapeContainsSubstroke(m_BestMatchingShape, stroke))
                {
                    string type = TopResult.GetShapePart(stroke);
                    ErrorDetail detail = ImageTemplateResult.GetDetail(type);
                    errors.Add(new ImageMatchError(ErrorType.Extra, detail, ErrorSeverity.Low, stroke));
                }
            }
            foreach (Substroke stroke in m_BestMatchingShape.SubstrokesL)
            {
                if (!Utilities.General.ShapeContainsSubstroke(m_GroupedShape, stroke))
                {
                    ErrorDetail detail = ImageTemplateResult.GetDetail(stroke.Labels[stroke.Labels.Length - 1]);
                    errors.Add(new ImageMatchError(ErrorType.Missing, detail, ErrorSeverity.Low, stroke));
                }
            }

            return errors;
        }

        public string ExpectedShapeName
        {
            get { return m_BestMatchingShape.Type; }
        }

        public bool ExpectedShapeIsParent
        {
            get { return (m_BestMatchingShapeParent == null); }
        }

        public string RecognizedShapeName
        {
            get 
            {
                if (TopResult != null)
                    return TopResult.Name;
                else
                    return "Unknown";
            }
        }

        public bool IsCorrectType
        {
            get
            {
                if (TopResult != null)
                    return (m_BestMatchingShape.Type == TopResult.Name);
                else
                    return false;
            }
        }

        public bool IsCorrectTypeInTop3
        {
            get
            {
                foreach (ImageTemplateResult result in m_ImageResults)
                {
                    if (result.Name == m_BestMatchingShape.Type)
                        return true;
                }

                return false;
            }
        }

        public bool IsCorrectCompleteness
        {
            get
            {
                if (TopResult == null)
                    return false;

                if (m_ExpectedErrors.Count != TopResult.Errors.Count)
                    return false;

                foreach (ImageMatchError expectedError in m_ExpectedErrors)
                {
                    bool foundMatch = false;
                    foreach (ImageMatchError error in TopResult.Errors)
                    {
                        if (IsEquivalentError(expectedError, error))
                        {
                            foundMatch = true;
                            continue;
                        }
                    }

                    if (!foundMatch)
                        return false;
                }

                return true;
            }
        }

        public bool DoesTopResultGivePerfectGuidance
        {
            get
            {
                if (TopResult == null)
                    return false;

                if (m_ExpectedErrors.Count != TopResult.Errors.Count)
                    return false;

                foreach (ImageMatchError expectedError in m_ExpectedErrors)
                {
                    bool foundMatch = false;
                    foreach (ImageMatchError error in TopResult.Errors)
                    {
                        if (IsSameError(expectedError, error))
                        {
                            foundMatch = true;
                            continue;
                        }
                    }

                    if (!foundMatch)
                        return false;
                }

                return true;
            }
        }

        private bool IsEquivalentError(ImageMatchError expectedError, ImageMatchError error)
        {
            if (expectedError.Type != error.Type)
                return false;

            return true;
        }

        private bool IsSameError(ImageMatchError expectedError, ImageMatchError error)
        {
            if (expectedError.Type != error.Type)
                return false;

            if (error.Type == ErrorType.Extra && expectedError.OffendingStroke.Id != error.OffendingStroke.Id)
                return false;

            if (error.Type == ErrorType.Missing && expectedError.Detail != error.Detail)
                return false;

            return true;
        }

        #region ISerializable Members

        public ShapeResult(SerializationInfo info, StreamingContext context)
        {
            m_BestMatchingShape = (Shape)info.GetValue("BestMatchingShape", typeof(Shape));
            m_BestMatchingShapeParent = (Shape)info.GetValue("BestMatchingShapeParent", typeof(Shape));
            m_ExpectedErrors = (List<ImageMatchError>)info.GetValue("ExpectedErrors", typeof(List<ImageMatchError>));
            m_GroupedShape = (Shape)info.GetValue("GroupedShape", typeof(Shape));
            m_ImageResults = (List<ImageTemplateResult>)info.GetValue("ImageResults", typeof(List<ImageTemplateResult>));
        }

        public void GetObjectData(SerializationInfo info, StreamingContext context)
        {
            info.AddValue("BestMatchingShape", m_BestMatchingShape);
            info.AddValue("BestMatchingShapeParent", m_BestMatchingShapeParent);
            info.AddValue("ExpectedErrors", m_ExpectedErrors);
            info.AddValue("GroupedShape", m_GroupedShape);
            info.AddValue("ImageResults", m_ImageResults);
        }

        #endregion

    }

    

    [Serializable]
    [DebuggerDisplay("Grouping: {m_GroupingCorrectness}, Completeness: {m_ResultCompleteness}, Confidence: {m_ResultConfidence}, Type: {m_ResultType}")]
    class ShapeCorrectness
    {
        GroupingCorrectness m_GroupingCorrectness;
        ResultCompleteness m_ResultCompleteness;
        ResultConfidence m_ResultConfidence;
        ResultType m_ResultType;

        public ShapeCorrectness()
        {
        }

        public GroupingCorrectness Grouping_Correctness
        {
            get { return m_GroupingCorrectness; }
            set { m_GroupingCorrectness = value; }
        }

        public ResultCompleteness Result_Completeness
        {
            get { return m_ResultCompleteness; }
            set { m_ResultCompleteness = value; }
        }

        public ResultConfidence Result_Confidence
        {
            get { return m_ResultConfidence; }
            set { m_ResultConfidence = value; }
        }

        public ResultType Result_Type
        {
            get { return m_ResultType; }
            set { m_ResultType = value; }
        }
    }
}
