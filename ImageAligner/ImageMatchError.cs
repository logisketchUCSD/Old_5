using System;
using System.Collections.Generic;
using System.Text;
using System.Diagnostics;

namespace ImageAligner
{
    public enum ErrorSeverity { High, Medium, Low };
    public enum ErrorType { Missing, Extra };
    public enum ErrorDetail {
        BackLine, BackArc, BottomArc, BottomLine, Bubble, FrontArc, 
        GreaterThan, Junk, Label, LabelBoxOther, Not_Hat, Not_V,
        TopArc, TopLine, TouchUp, Triangle, Unknown, Wire };

    [DebuggerDisplay("{m_Type} {m_Detail} - {m_Severity}")]
    [Serializable]
    public class ImageMatchError
    {
        /// <summary>
        /// i.e. 'Missing' or 'Extra'
        /// </summary>
        ErrorType m_Type;

        /// <summary>
        /// e.g. BackArc (modifying 'Missing' type)
        /// </summary>
        ErrorDetail m_Detail;

        /// <summary>
        /// 'High', 'Medium', 'Low'
        /// </summary>
        ErrorSeverity m_Severity;

        Sketch.Substroke m_SubstrokeForError;

        public ImageMatchError(ErrorType type, ErrorDetail detail, ErrorSeverity severity)
        {
            m_Type = type;
            m_Detail = detail;
            m_Severity = severity;
            m_SubstrokeForError = null;
        }

        public ImageMatchError(ErrorType type, ErrorDetail detail, ErrorSeverity severity, Sketch.Substroke offendingStroke)
        {
            m_Type = type;
            m_Detail = detail;
            m_Severity = severity;
            m_SubstrokeForError = offendingStroke;
        }

        public ErrorType Type
        {
            get { return m_Type; }
        }

        public ErrorDetail Detail
        {
            get { return m_Detail; }
        }

        public ErrorSeverity Severity
        {
            get { return m_Severity; }
        }

        public Sketch.Substroke OffendingStroke
        {
            get { return m_SubstrokeForError; }
        }
    }
}
