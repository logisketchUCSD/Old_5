/*
 * File: SymbolInfo.cs
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
using System.Diagnostics;
using Utilities;

namespace ImageAligner
{
    /// <summary>
    /// This object simply holds general information about a symbol: the 
    /// type (name), class of the symbol, user who drew it, completeness, 
    /// platform drawn on, task when drawing, and a unique id.
    /// </summary>
    [Serializable]
    [DebuggerDisplay("Type={m_SymbolType}, User={m_User.Name}")]
    public class SymbolInfo
    {
        #region Member Variables

        /// <summary>
        /// Unique ID for this symbol
        /// </summary>
        Guid m_SymbolId;

        /// <summary>
        /// The type of symbol this ImageTemplate describes (e.g. AND, OR, Arrow, etc.)
        /// </summary>
        string m_SymbolType;

        /// <summary>
        /// The symbol's container class (e.g. AND --> Gate, Wire --> Connector)
        /// </summary>
        string m_SymbolClass;

        /// <summary>
        /// The user who drew this shape
        /// </summary>
        User m_User;

        /// <summary>
        /// The platform this shape was drawn on (i.e. TabletPC or Wacom Tablet (paper))
        /// </summary>
        PlatformUsed m_PlatformUsed;

        /// <summary>
        /// How complete this symbol is
        /// </summary>
        SymbolCompleteness m_Completeness;

        /// <summary>
        /// Under what circumstances was this shape drawn
        /// </summary>
        DrawingTask m_DrawingTask;

        #endregion

        #region Constructors

        /// <summary>
        /// Default Constructor, sets values to defaults or "None"
        /// </summary>
        public SymbolInfo()
        {
            m_SymbolId = Guid.NewGuid();
            m_SymbolType = "None";
            m_SymbolClass = "None";
            m_User = new User();
            m_PlatformUsed = PlatformUsed.TabletPC;
            m_Completeness = SymbolCompleteness.Complete;
            m_DrawingTask = DrawingTask.Synthesize;
        }

        /// <summary>
        /// Constructor with limited information
        /// </summary>
        /// <param name="name"></param>
        public SymbolInfo(string name)
        {
            m_SymbolId = Guid.NewGuid();
            m_SymbolType = name;
            m_SymbolClass = "Unknown";
            m_User = new User();
            m_PlatformUsed = PlatformUsed.TabletPC;
            m_Completeness = SymbolCompleteness.Complete;
            m_DrawingTask = DrawingTask.Synthesize;
        }

        /// <summary>
        /// Typical Constructor
        /// </summary>
        /// <param name="user">User who drew the symbol</param>
        /// <param name="symbolType">Type of symbol drawn</param>
        /// <param name="symbolClass">Class of drawn symbol</param>
        public SymbolInfo(User user, string symbolType, string symbolClass)
        {
            m_SymbolId = Guid.NewGuid();
            m_SymbolType = symbolType;
            m_SymbolClass = symbolClass;
            m_User = user;
            m_PlatformUsed = PlatformUsed.TabletPC;
            m_Completeness = SymbolCompleteness.Complete;
            m_DrawingTask = DrawingTask.Synthesize;
        }

        #endregion

        #region Getters/Setters

        /// <summary>
        /// Unique ID for this symbol
        /// </summary>
        public Guid SymbolId
        {
            get { return m_SymbolId; }
            set { m_SymbolId = value; }
        }

        /// <summary>
        /// The type of symbol this ImageTemplate describes (e.g. AND, OR, Arrow, etc.)
        /// </summary>
        public string SymbolType
        {
            get { return m_SymbolType; }
            set { m_SymbolType = value; }
        }

        /// <summary>
        /// The symbol's container class (e.g. AND --> Gate, Wire --> Connector)
        /// </summary>
        public string SymbolClass
        {
            get { return m_SymbolClass; }
            set { m_SymbolClass = value; }
        }

        /// <summary>
        /// The user who drew this shape
        /// </summary>
        public User User
        {
            get { return m_User; }
            set { m_User = value; }
        }

        /// <summary>
        /// The platform this shape was drawn on (i.e. TabletPC or Wacom Tablet (paper))
        /// </summary>
        public PlatformUsed PlatformUsed
        {
            get { return m_PlatformUsed; }
            set { m_PlatformUsed = value; }
        }

        /// <summary>
        /// How complete this symbol is
        /// </summary>
        public SymbolCompleteness Completeness
        {
            get { return m_Completeness; }
            set { m_Completeness = value; }
        }

        /// <summary>
        /// Under what circumstances was this shape drawn
        /// </summary>
        public DrawingTask DrawTask
        {
            get { return m_DrawingTask; }
            set { m_DrawingTask = value; }
        }

        #endregion
    }
}
