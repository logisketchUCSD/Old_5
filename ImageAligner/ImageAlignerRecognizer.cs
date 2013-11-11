using System;
using System.Collections.Generic;
using System.Text;
using System.Runtime.Serialization;
using System.Runtime.Serialization.Formatters.Binary;
using Sketch;
using Utilities;

namespace ImageAligner
{
    [Serializable]
    public class ImageAlignerRecognizer : ISerializable
    {
        #region Member Variables

        /// <summary>
        /// Single-Stroke Recognizer
        /// </summary>
        RubineRecognizerUpdateable m_Rubine;

        bool m_NeedsUpdating;

        /// <summary>
        /// ImageTemplate List
        /// </summary>
        List<ImageTemplate> m_Templates;

        User m_User;

        PlatformUsed m_Platform;

        ConfusionMatrix m_SSRConfusionMatrix;
        string c_ConfusionMatrixFile = "Code\\Recognition\\ImageAligner\\allConfusion.txt";

        #endregion

        #region Constructors

        public ImageAlignerRecognizer()
        {
            m_NeedsUpdating = false;
            m_Rubine = new RubineRecognizerUpdateable();
            m_Templates = new List<ImageTemplate>();
            m_User = new User();
            m_Platform = PlatformUsed.TabletPC;
            GetConfusionMatrix(c_ConfusionMatrixFile);
        }

        public ImageAlignerRecognizer(User user, PlatformUsed platform)
        {
            m_NeedsUpdating = false;
            m_Rubine = new RubineRecognizerUpdateable();
            m_Templates = new List<ImageTemplate>();
            m_User = user;
            m_Platform = platform;
            GetConfusionMatrix(c_ConfusionMatrixFile);
        }

        private void GetConfusionMatrix(string filename)
        {
            string filePath = System.IO.Directory.GetCurrentDirectory();
            if (filePath.Contains("\\Code\\"))
                filePath = filePath.Substring(0, filePath.IndexOf("\\Code\\") + 1);
            else if (filePath.Contains("\\Sketch\\"))
                filePath = filePath.Substring(0, filePath.IndexOf("\\Sketch\\") + 8);
            else if (filePath.Contains("\\sketch\\"))
                filePath = filePath.Substring(0, filePath.IndexOf("\\sketch\\") + 8);
            else if (filePath.Contains("\\Trunk\\"))
                filePath = filePath.Substring(0, filePath.IndexOf("\\Trunk\\") + 7);
            filePath += c_ConfusionMatrixFile;
            m_SSRConfusionMatrix = new ConfusionMatrix(General.leafLabels);
            m_SSRConfusionMatrix.LoadFromFile(filePath);
        }

        #endregion

        #region Recognition

        public List<ImageTemplateResult> Recognize(Shape shape, int n)
        {
            if (m_NeedsUpdating)
            {
                bool updated = m_Rubine.updateMatrices();
                if (updated)
                    m_NeedsUpdating = false;
            }

            if (m_NeedsUpdating)
                return null;

            Dictionary<Substroke, string> rubineResults = RecognizeSingleStrokes(shape);

            ImageTemplate unknown = new ImageTemplate(shape, rubineResults, m_SSRConfusionMatrix);

            List<ImageTemplateResult> results = unknown.Recognize(m_Templates, n);

            return results;
        }

        public ImageTemplateResult Recognize(Shape shape)
        {
            if (m_NeedsUpdating)
            {
                bool updated = m_Rubine.updateMatrices();
                if (updated)
                    m_NeedsUpdating = false;
            }

            if (m_NeedsUpdating)
                return null;

            Dictionary<Substroke, string> rubineResults = RecognizeSingleStrokes(shape);

            ImageTemplate unknown = new ImageTemplate(shape, rubineResults, m_SSRConfusionMatrix);

            ImageTemplateResult result = unknown.Recognize(m_Templates);

            return result;
        }

        private Dictionary<Substroke, string> RecognizeSingleStrokes(Shape shape)
        {
            Dictionary<Substroke, string> parts = new Dictionary<Substroke, string>();
            foreach (Substroke stroke in shape.Substrokes)
                parts.Add(stroke, m_Rubine.Recognize(stroke));

            return parts;
        }

        #endregion

        #region Other Functions

        public void Add(Shape shape)
        {
            Dictionary<Substroke, string> parts = new Dictionary<Substroke, string>();
            foreach (Substroke s in shape.Substrokes)
                parts.Add(s, s.Labels[s.Labels.Length - 1]);

            SymbolInfo symbolInfo = new SymbolInfo();
            symbolInfo.Completeness = SymbolCompleteness.Complete;
            symbolInfo.SymbolClass = General.GetClass(shape);
            symbolInfo.SymbolType = shape.Type;
            symbolInfo.User = m_User;
            symbolInfo.PlatformUsed = m_Platform;

            ImageTemplate template = new ImageTemplate(shape, parts, symbolInfo, m_SSRConfusionMatrix);
            m_Templates.Add(template);

            foreach (Substroke s in shape.Substrokes)
                m_Rubine.Add(s.Labels[s.Labels.Length - 1], s.PointsL);

            m_NeedsUpdating = true;
        }

        #endregion

        
        #region Serialization

        /// <summary>
        /// Deserialization Constructor
        /// </summary>
        /// <param name="info"></param>
        /// <param name="ctxt"></param>
        public ImageAlignerRecognizer(SerializationInfo info, StreamingContext ctxt)
        {
            //Get the values from info and assign them to the appropriate properties
            m_Rubine = (RubineRecognizerUpdateable)info.GetValue("Rubine", typeof(RubineRecognizerUpdateable));
            m_Templates = (List<ImageTemplate>)info.GetValue("Templates", typeof(List<ImageTemplate>));
            m_User = (User)info.GetValue("User", typeof(User));
            m_Platform = (PlatformUsed)info.GetValue("Platform", typeof(PlatformUsed));
            bool foundMatrix = false;
            try
            {
                m_SSRConfusionMatrix = (ConfusionMatrix)info.GetValue("ConfusionMatrix", typeof(ConfusionMatrix));
                foundMatrix = true;
            }
            catch (Exception e)
            {
                Console.WriteLine(e.Message);
            }
            if (!foundMatrix)
                GetConfusionMatrix(c_ConfusionMatrixFile);

            m_NeedsUpdating = true;
        }

        /// <summary>
        /// Serialization Function
        /// </summary>
        /// <param name="info"></param>
        /// <param name="ctxt"></param>
        public void GetObjectData(SerializationInfo info, StreamingContext ctxt)
        {
            info.AddValue("Rubine", m_Rubine);
            info.AddValue("Templates", m_Templates);
            info.AddValue("User", m_User);
            info.AddValue("Platform", m_Platform);
            info.AddValue("ConfusionMatrix", m_SSRConfusionMatrix);
        }

        #endregion

        public static ImageAlignerRecognizer Load(string filename)
        {
            System.IO.Stream stream = System.IO.File.Open(filename, System.IO.FileMode.Open);
            BinaryFormatter bformatter = new BinaryFormatter();
            ImageAlignerRecognizer aligner = (ImageAlignerRecognizer)bformatter.Deserialize(stream);
            stream.Close();

            return aligner;
        }

        public void Save(string filename)
        {
            System.IO.Stream stream = System.IO.File.Open(filename, System.IO.FileMode.Create);
            BinaryFormatter bformatter = new BinaryFormatter();
            bformatter.Serialize(stream, this);
            stream.Close();
        }
    }
}
