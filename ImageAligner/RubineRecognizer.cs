using System;
using System.IO;
using System.Collections.Generic;
using System.Text;
using System.Runtime.Serialization;
using Utilities.Matrix;
using Sketch;

namespace ImageAligner
{
    /// <summary>
    /// This is a lightweight version of the Rubine Recognizer 
    /// which can only classify/recognize new instances, it cannot 
    /// be trained.
    /// </summary>
    [Serializable]
    public class RubineRecognizer : ISerializable
    {
        private const double SCORE_THRESHOLD = 0.95;

        #region Member Variables
        GeneralMatrix _weights;
        List<double> _weights_0;
        List<string> _classNames;
        #endregion

        #region Constructors
        public RubineRecognizer(GeneralMatrix weights, List<double> weights_0, List<string> classNames)
        {
            _weights = weights;
            _weights_0 = weights_0;
            _classNames = classNames;
        }

        public RubineRecognizer()
        {
            _weights = new GeneralMatrix(1, 1);
            _weights_0 = new List<double>();
            _classNames = new List<string>();
        }
        #endregion

        #region Interface Functions
        public string Recognize(Sketch.Substroke stroke)
        {
            List<double> features = getFeatureValues(stroke.PointsL);
            string best = "NA";
            double bestScore = double.MinValue;
            List<double> allScores = new List<double>();

            // Compute scores for each class
            for (int c = 0; c < _classNames.Count; c++)
            {
                double sum = 0.0;
                for (int i = 0; i < _weights_0.Count; i++)
                    sum += _weights.GetElement(c, i) * features[i];

                double score = this._weights_0[c] + sum;
                allScores.Add(score);
                if (score > bestScore)
                {
                    bestScore = score;
                    best = _classNames[c];
                }
            }
            allScores.Remove(bestScore);

            // Get Probability - for rejection
            double sumP = 0.0;
            foreach (double s in allScores)
                    sumP += Math.Exp(s - bestScore);
            double probability = 1 / sumP;

            if (probability < SCORE_THRESHOLD)
                return "Unknown";
            //stroke.StrokeType = best;
            return best;
        }

        private List<double> getFeatureValues(List<Point> points)
        {
            RubineTemplate temp = new RubineTemplate(points);
            return temp.Features;
        }

        #endregion

        #region Serialization

        /// <summary>
        /// Deserialization Constructor
        /// </summary>
        /// <param name="info"></param>
        /// <param name="ctxt"></param>
        public RubineRecognizer(SerializationInfo info, StreamingContext ctxt)
        {
            //Get the values from info and assign them to the appropriate properties
            _weights = (GeneralMatrix)info.GetValue("weights", typeof(GeneralMatrix));
            _weights_0 = (List<double>)info.GetValue("weights_0", typeof(List<double>));
            _classNames = (List<string>)info.GetValue("classes", typeof(List<string>));
        }

        /// <summary>
        /// Serialization Function
        /// </summary>
        /// <param name="info"></param>
        /// <param name="ctxt"></param>
        public void GetObjectData(SerializationInfo info, StreamingContext ctxt)
        {
            info.AddValue("weights", _weights);
            info.AddValue("weights_0", _weights_0);
            info.AddValue("classes", _classNames);
        }

        #endregion


        public static RubineRecognizer Load(string filename)
        {
            System.IO.Stream stream = System.IO.File.Open(filename, System.IO.FileMode.Open);
            System.Runtime.Serialization.Formatters.Binary.BinaryFormatter bformatter = new System.Runtime.Serialization.Formatters.Binary.BinaryFormatter();
            RubineRecognizer rubine = (RubineRecognizer)bformatter.Deserialize(stream);
            stream.Close();

            return rubine;
        }

        public void Save(string filename)
        {
            System.IO.Stream stream = System.IO.File.Open(filename, System.IO.FileMode.Create);
            System.Runtime.Serialization.Formatters.Binary.BinaryFormatter bformatter = new System.Runtime.Serialization.Formatters.Binary.BinaryFormatter();
            bformatter.Serialize(stream, this);
            stream.Close();
        }
    }

    /// <summary>
    /// This Rubine Recognizer has all the feature values for 
    /// every example stored in it, making it quite large. A 
    /// RubineRecognizer can be created if a lightweight version 
    /// is needed, assuming you have already trained this one.
    /// </summary>
    [Serializable]
    public class RubineRecognizerUpdateable : ISerializable
    {
        private const double SCORE_THRESHOLD = 0.95;

        #region Member Variables
        GeneralMatrix ccm;
        GeneralMatrix ccmInverse;
        GeneralMatrix weights;
        List<double> weights_0;
        List<RubineClass> classes;
        #endregion

        #region Constructors
        public RubineRecognizerUpdateable()
        {
            ccm = new GeneralMatrix(1, 1);
            ccmInverse = new GeneralMatrix(1, 1);
            weights = new GeneralMatrix(1, 1);
            weights_0 = new List<double>();
            classes = new List<RubineClass>();
        }
        #endregion

        #region Getters & Setters
        public GeneralMatrix CCM
        { 
            get { return this.ccm; }
        }

        public GeneralMatrix CCMinverse
        {
            get { return this.ccmInverse; }
        }

        public GeneralMatrix Weights
        {
            get { return this.weights; }
        }

        public List<double> Weights_0
        {
            get { return this.weights_0; }
        }

        public List<RubineClass> Classes
        {
            get { return this.classes; }
        }

        public int NumExamples
        {
            get 
            {
                int numExamples = 0;
                for (int i = 0; i < this.classes.Count; i++)
                {
                    numExamples += classes[i].Examples.Count;
                }
                return numExamples;
            }
        }

        public RubineRecognizer LiteRecognizer
        {
            get
            {
                List<string> cls = new List<string>(classes.Count);
                foreach (RubineClass c in classes)
                    cls.Add(c.Name);

                RubineRecognizer r = new RubineRecognizer(weights, weights_0, cls);

                return r;
            }
        }

        #endregion

        #region Interface Functions
        public void Add(string name, List<Sketch.Point> points)
        {
            List<string> allNames = new List<string>();
            foreach (RubineClass c in classes)
                allNames.Add(c.Name);

            if (allNames.Contains(name))
                addToClass(name, points);
            else
                addAsNewClass(name, points);
        }

        public void addAsNewClass(string name, List<Sketch.Point> points)
        {
            classes.Add(new RubineClass(name, points));
            if (classes.Count == 1)
                updateMatrixSize();
        }

        public void addToClass(string name, List<Sketch.Point> points)
        {
            for (int i = 0; i < classes.Count; i++)
            {
                if (string.Equals(name, classes[i].Name))
                    classes[i].addExample(points);
            }
            //calculateCCM();
        }

        /*public string loadTrainingData()
        {
            //string dir = "C:/Documents and Settings/eric/My Documents/CURRENT/Classes - ME231 - Pen-Based Computing/Project 3 - Rubine Recognizer/Training Data/";
            OpenFileDialog dlg = new OpenFileDialog();
            dlg.Filter = "Text Files (*.txt)|*.txt";
            dlg.Multiselect = true;
            dlg.RestoreDirectory = true;
            //dlg.InitialDirectory = "C:\\Documents and Settings\\eric\\My Documents\\CURRENT\\Classes - ME231 - Pen-Based Computing\\Project 3 - Rubine Recognizer\\Training Data\\";

            string directory = "C:";
            if (dlg.ShowDialog() == DialogResult.OK)
            {
                directory = Path.GetDirectoryName(dlg.FileNames[0]);
                for (int i = 0; i < dlg.FileNames.Length; i++)
                    readXY(dlg.FileNames[i]);
            }

            return directory;
        }*/

        public bool updateMatrices()
        {
            if (classes.Count > 1)
            {
                foreach (RubineClass c in classes)
                {
                    //if (c.Examples.Count <= 1)
                        //return false;
                    c.calculateAverages();
                    c.calculateCovarianceMatrix();
                }

                calculateCCM();
                calculateCCMInverse();
                updateWeightsSize();
                calculateWeights();
                calculateWeights_0();

                return true;
            }
            return false;
        }

        public string Recognize(Sketch.Substroke stroke)
        {
            List<double> features = getFeatureValues(stroke.PointsL);
            string best = "NA";
            double bestScore = double.MinValue;
            List<double> allScores = new List<double>();

            // Compute scores for each class
            for (int c = 0; c < classes.Count; c++)
            {
                double sum = 0.0;
                for (int i = 0; i < classes[c].Averages.Count; i++)
                    sum += weights.GetElement(c, i) * features[i];

                double score = this.weights_0[c] + sum;
                allScores.Add(score);
                if (score > bestScore)
                {
                    bestScore = score;
                    best = classes[c].Name;
                }
            }
            allScores.Remove(bestScore);

            // Get Probability - for rejection
            double sumP = 0.0;
            foreach (double s in allScores)
                    sumP += Math.Exp(s - bestScore);
            double probability = 1 / sumP;

            //if (probability < SCORE_THRESHOLD)
                //return "Unknown";

            return best;
        }

        private List<double> getFeatureValues(List<Point> points)
        {
            RubineTemplate temp = new RubineTemplate(points);
            return temp.Features;
        }

        public List<String> classifyGesture(List<Sketch.Point> points, ref double probability)
        {
            List<String> classifiedName = new List<String>();
            double sumI = new double();
            List<double> features = this.classes[0].Examples[0].getFeatures(points);

            List<double> vec = new List<double>(this.classes.Count);

            for (int c = 0; c < this.classes.Count; c++)
            {
                sumI = 0.0;
                for (int i = 0; i < classes[c].Averages.Count; i++)
                    sumI += weights.GetElement(c, i) * features[i];

                double value = this.weights_0[c] + sumI;
                vec.Add(value);
            }

            double best = double.MinValue;
            int loc = 0;
            for (int i = 0; i < vec.Count; i++)
            {
                if (vec[i] > best)
                {
                    best = vec[i];
                    loc = i;
                    //classifiedName.Add(this.classes[i].Name);
                }
            }

            double sumP = 0.0;
            for (int i = 0; i < vec.Count; i++)
            {
                if (i != loc)
                {
                    sumP += Math.Exp((vec[i] - vec[loc]) / vec[loc]);
                }
            }
            probability = 1 / sumP;

            List<double> newVec = new List<double>(vec);
            newVec.Sort();

            for (int i = newVec.Count - 1; i >= 0; i--)
            {
                for (int j = 0; j < vec.Count; j++)
                {
                    if (vec[j] == newVec[i])
                        classifiedName.Add(this.classes[j].Name);
                }
            }

            return classifiedName;
        }
        #endregion

        #region Private Functions
        private void updateMatrixSize()
        {
            this.ccm = new GeneralMatrix(this.classes[0].CovarianceMatrix.RowDimension, this.classes[0].CovarianceMatrix.ColumnDimension);
            this.ccmInverse = new GeneralMatrix(this.classes[0].CovarianceMatrix.RowDimension, this.classes[0].CovarianceMatrix.ColumnDimension);
            this.weights = new GeneralMatrix(this.classes.Count, this.classes[0].CovarianceMatrix.RowDimension);
            this.weights_0 = new List<double>(this.classes.Count);
        }

        private void updateWeightsSize()
        {
            this.weights = new GeneralMatrix(this.classes.Count, this.classes[0].CovarianceMatrix.RowDimension, 0.0);
            
            this.weights_0 = new List<double>(this.classes.Count);
        }

        private void calculateCCM()
        {
            double sumC = new double();
            int numClasses = this.classes.Count;
            int numExamples = 0;
            for (int i = 0; i < this.classes.Count; i++)
                numExamples += this.classes[i].Examples.Count;
            int denom = numExamples - numClasses;


            for (int i = 0; i < classes[0].Averages.Count; i++)
            {
                for (int j = 0; j < classes[0].Averages.Count; j++)
                {
                    sumC = 0.0;
                    for (int c = 0; c < classes.Count; c++)
                        sumC += classes[c].CovarianceMatrix.GetElement(i, j);

                    if (denom > 0)
                        this.ccm.SetElement(i, j, sumC / (double)denom);
                    else
                        this.ccm.SetElement(i, j, 0.0);
                }
            }
        }

        private void calculateCCMInverse()
        {
            this.ccmInverse = ccm.Inverse();
        }

        private void calculateWeights()
        {
            double sumI = new double();
            for (int c = 0; c < this.classes.Count; c++)
            {
                for (int j = 0; j < this.classes[c].Averages.Count; j++)
                {
                    sumI = 0.0;
                    for (int i = 0; i < this.classes[c].Averages.Count; i++)
                    {
                        sumI += this.ccmInverse.GetElement(i, j) * this.classes[c].Averages[i];
                    }
                    this.weights.SetElement(c, j, sumI);
                }
            }
        }

        private void calculateWeights_0()
        {
            double sumI = new double();

            for (int c = 0; c < classes.Count; c++)
            {
                sumI = 0.0;
                for (int i = 0; i < this.classes[c].Averages.Count; i++)
                {
                    sumI += this.weights.GetElement(c, i) * this.classes[c].Averages[i];
                }

                this.weights_0.Add(-0.5 * sumI);
            }
        }

        private void readXY(string filename)
        {
            string name = "";
            string strLine = "";
            string[] strs = new string[0];
            int numPoints = 0;
            Sketch.Point q = new Sketch.Point();

            StreamReader reader = new StreamReader(filename);

            name = reader.ReadLine();
            strLine = reader.ReadLine();
            numPoints = Convert.ToInt32(strLine);
            List<Sketch.Point> points = new List<Sketch.Point>(numPoints);
            float X = 0;
            float Y = 0;
            ulong T = 0;

            for (int i = 0; i < numPoints; i++)
            {
                strLine = reader.ReadLine();
                strs = strLine.Split(new char[] { ',' } , 3);
                X = (float)Convert.ToInt32(strs[0]);
                Y = (float)Convert.ToInt32(strs[1]);
                T = Convert.ToUInt64(strs[2]);
                Sketch.Point p = new Sketch.Point(X, Y);
                p.Time = T;
                points.Add(p);
            }

            bool newClass = true;
            for (int i = 0; i < this.classes.Count; i++)
            {
                if (name == this.classes[i].Name)
                    newClass = false;
            }

            if (newClass)
                addAsNewClass(name, points);
            else
                addToClass(name, points);

            reader.Close();
        }
        #endregion

        #region Serialization

        /// <summary>
        /// Deserialization Constructor
        /// </summary>
        /// <param name="info"></param>
        /// <param name="ctxt"></param>
        public RubineRecognizerUpdateable(SerializationInfo info, StreamingContext ctxt)
        {
            //Get the values from info and assign them to the appropriate properties
            ccm = (GeneralMatrix)info.GetValue("ccm", typeof(GeneralMatrix));
            ccmInverse = (GeneralMatrix)info.GetValue("ccmInverse", typeof(GeneralMatrix));
            weights = (GeneralMatrix)info.GetValue("weights", typeof(GeneralMatrix));
            weights_0 = (List<double>)info.GetValue("weights_0", typeof(List<double>));
            classes = (List<RubineClass>)info.GetValue("classes", typeof(List<RubineClass>));
        }

        /// <summary>
        /// Serialization Function
        /// </summary>
        /// <param name="info"></param>
        /// <param name="ctxt"></param>
        public void GetObjectData(SerializationInfo info, StreamingContext ctxt)
        {
            info.AddValue("ccm", ccm);
            info.AddValue("ccmInverse", ccmInverse);
            info.AddValue("weights", weights);
            info.AddValue("weights_0", weights_0);
            info.AddValue("classes", classes);
        }

        #endregion

        public static RubineRecognizerUpdateable Load(string filename)
        {
            System.IO.Stream stream = System.IO.File.Open(filename, System.IO.FileMode.Open);
            System.Runtime.Serialization.Formatters.Binary.BinaryFormatter bformatter = new System.Runtime.Serialization.Formatters.Binary.BinaryFormatter();
            RubineRecognizerUpdateable rubine = (RubineRecognizerUpdateable)bformatter.Deserialize(stream);
            stream.Close();

            return rubine;
        }

        public void Save(string filename)
        {
            System.IO.Stream stream = System.IO.File.Open(filename, System.IO.FileMode.Create);
            System.Runtime.Serialization.Formatters.Binary.BinaryFormatter bformatter = new System.Runtime.Serialization.Formatters.Binary.BinaryFormatter();
            bformatter.Serialize(stream, this);
            stream.Close();
        }
    }

    [Serializable]
    public class RubineClass
    {
        #region Member Variables
        string name;
        List<RubineTemplate> examples;
        List<double> averages;
        GeneralMatrix covarianceMatrix;
        #endregion

        #region Constructors
        public RubineClass()
        {
            name = "NA";
            examples = new List<RubineTemplate>();
            averages = new List<double>();
            covarianceMatrix = new GeneralMatrix(1, 1);
        }

        public RubineClass(string name, List<Sketch.Point> points)
        {
            this.name = name;
            this.examples = new List<RubineTemplate>();
            this.examples.Add(new RubineTemplate(points));
            calculateAverages();
            this.covarianceMatrix = new GeneralMatrix(this.averages.Count, this.averages.Count, 0.0);
        }
        #endregion

        #region Getters & Setters
        public string Name
        {
            get { return this.name; }
            set { this.name = value; }
        }

        public List<RubineTemplate> Examples
        {
            get { return this.examples; }
        }

        public GeneralMatrix CovarianceMatrix
        {
            get { return this.covarianceMatrix; }
        }

        public List<double> Averages
        {
            get { return this.averages; }
        }
        #endregion

        #region Interface Functions
        public void addExample(List<Sketch.Point> points)
        {
            this.examples.Add(new RubineTemplate(points));
            //calculateAverages();
            //calculateCovarianceMatrix();
        }
        #endregion

        #region Private Functions
        public void calculateAverages()
        {
            this.averages = new List<double>();
            
            for (int i = 0; i < this.examples[0].Features.Count; i++)
            {
                double sumFcei = 0.0;

                for (int e = 0; e < this.examples.Count; e++)
                    sumFcei += this.examples[e].Features[i];

                this.averages.Add(sumFcei / this.examples.Count);
            }
        }

        /// <summary>
        /// Calculates the covariance matrix within a given class
        /// </summary>
        public void calculateCovarianceMatrix()
        {
            if (this.covarianceMatrix == null || this.covarianceMatrix.RowDimension != this.averages.Count)
                this.covarianceMatrix = new GeneralMatrix(this.averages.Count, this.averages.Count, 0.0);

            int numFeatures = this.examples[0].Features.Count;
            int numExamples = this.examples.Count;
            for (int i = 0; i < numFeatures; i++)
            {
                for (int j = 0; j < numFeatures; j++)
                {
                    double sumE = 0.0;
                    for (int e = 0; e < numExamples; e++)
                    {
                        double diff_i = this.examples[e].Features[i] - this.averages[i];
                        double diff_j = this.examples[e].Features[j] - this.averages[j];
                        sumE += diff_i * diff_j;
                    }
                    this.covarianceMatrix.SetElement(i, j, sumE);
                }
            }
        }
        #endregion
    }

    [Serializable]
    public class RubineTemplate
    {
        #region Member Variables
        List<Sketch.Point> points;
        List<double> features;
        #endregion

        #region Constructors
        public RubineTemplate()
        {
            features = new List<double>();
            points = new List<Sketch.Point>();
        }

        public RubineTemplate(List<Sketch.Point> points)
        {
            this.points = points;
            this.features = getFeatures(points);
        }
        #endregion

        #region Getters & Setters
        public List<Sketch.Point> Points
        {
            get { return this.points; }
        }

        public List<double> Features
        {
            get { return this.features; }
        }
        #endregion

        #region Interface Functions
        public void makeRubineTemplate(List<Sketch.Point> points)
        {
            this.points = points;
            this.features = getFeatures(points);
        }

        public List<double> getFeatures(List<Sketch.Point> points)
        {
            List<double> features = new List<double>();
            System.Drawing.RectangleF bbox;

            features.Add(feature1(points));
            features.Add(feature2(points));

            features.Add(feature3(points, out bbox));
            features.Add(feature4(points, bbox));

            double f5 = feature5(points);
            features.Add(f5);
            features.Add(feature6(points, f5));
            features.Add(feature7(points, f5));

            features.Add(feature8(points));

            List<double> thetas = GetThetas(points);
            features.Add(feature9(thetas));
            features.Add(feature10(thetas));
            features.Add(feature11(thetas));

            features.Add(feature12(points));
            features.Add(feature13(points));

            return features;
        }

        #endregion

        #region Private Functions

        private List<double> GetThetas(List<Point> points)
        {
            List<double> thetas = new List<double>(points.Count);

            for (int i = 1; i < points.Count - 1; i++)
            {
                float deltax = points[i + 1].X - points[i].X;
                float deltay = points[i + 1].Y - points[i].Y;
                float deltax1 = points[i].X - points[i - 1].X;
                float deltay1 = points[i].Y - points[i - 1].Y;

                thetas.Add( Math.Atan2(
                    deltax * deltay1 - deltax1 * deltay,
                    deltax * deltax1 + deltay * deltay1));
            }

            return thetas;
        }

        private double feature1(List<Sketch.Point> points)
        {
            double featureValue = 0.0;

            if (points.Count > 2)
            {
                Sketch.Point p1 = points[0];
                int n = 2;
                Sketch.Point p2 = points[n];
                while (p1.X == p2.X && p1.Y == p2.Y)
                {
                    n++;
                    if (n < points.Count)
                        p2 = points[n];
                }

                double numer = (double)p2.X - (double)p1.X;
                double denom = Math.Sqrt(Math.Pow(p2.X - p1.X, 2.0) + Math.Pow(p2.Y - p1.Y, 2.0));
                
                if (denom != 0)
                    featureValue = numer / denom;
                else
                    featureValue = 0.0;
            }

            return featureValue;
        }

        private double feature2(List<Sketch.Point> points)
        {
            double featureValue = 0.0;

            if (points.Count > 2)
            {
                Sketch.Point p1 = points[0];
                int n = 2;
                Sketch.Point p2 = points[n];
                while (p1.X == p2.X && p1.Y == p2.Y)
                {
                    n++;
                    if (n < points.Count)
                        p2 = points[n];
                }

                double numer = (double)p2.Y - (double)p1.Y;
                double denom = Math.Sqrt(Math.Pow(p2.X - p1.X, 2.0) + Math.Pow(p2.Y - p1.Y, 2.0));

                if (denom != 0)
                    featureValue = numer / denom;
                else
                    featureValue = 0.0;
            }

            return featureValue;
        }

        private double feature3(List<Sketch.Point> points, out System.Drawing.RectangleF bbox)
        {
            double featureValue = 0.0;

            float maxX = float.MinValue;
            float minX = float.MaxValue;
            float maxY = float.MinValue;
            float minY = float.MaxValue;

            for (int i = 0; i < points.Count; i++)
            {
                if (points[i].X > maxX) maxX = points[i].X;
                if (points[i].X < minX) minX = points[i].X;
                if (points[i].Y > maxY) maxY = points[i].Y;
                if (points[i].Y < minY) minY = points[i].Y;
            }

            bbox = new System.Drawing.RectangleF(minX, minY, maxX - minX, maxY - minY);

            featureValue = Math.Sqrt(Math.Pow(maxX - minX, 2.0) + Math.Pow(maxY - minY, 2.0));

            return featureValue;
        }

        private double feature4(List<Sketch.Point> points, System.Drawing.RectangleF bbox)
        {
            double featureValue = 0.0;
            float maxX = bbox.X + bbox.Width;
            float minX = bbox.X;
            float maxY = bbox.Y + bbox.Height;
            float minY = bbox.Y;

            featureValue = Math.Atan2((double)maxY - (double)minY, (double)maxX - (double)minX);

            return featureValue;
        }

        private double feature5(List<Sketch.Point> points)
        {
            double featureValue = 0.0;

            featureValue = Math.Sqrt(
                Math.Pow(points[points.Count - 1].X - points[0].X, 2.0) +
                Math.Pow(points[points.Count - 1].Y - points[0].Y, 2.0) );

            return featureValue;
        }

        private double feature6(List<Sketch.Point> points, double f5)
        {
            double featureValue = 0.0;

            if (f5 != 0.0)
                featureValue = (points[points.Count - 1].X - points[0].X) / f5;

            return featureValue;
        }

        private double feature7(List<Sketch.Point> points, double f5)
        {
            double featureValue = 0.0;

            if (f5 != 0.0)
                featureValue = (points[points.Count - 1].Y - points[0].Y) / f5;

            return featureValue;
        }

        private double feature8(List<Sketch.Point> points)
        {
            double featureValue = 0.0;
            float deltax, deltay;

            for (int i = 0; i < points.Count - 1; i++)
            {
                deltax = points[i + 1].X - points[i].X;
                deltay = points[i + 1].Y - points[i].Y;
                featureValue += Math.Sqrt(
                    Math.Pow(deltax, 2.0) +
                    Math.Pow(deltay, 2.0));
            }

            return featureValue;
        }

        private double feature9(List<double> thetas)
        {
            double featureValue = 0.0;

            foreach (double t in thetas)
                featureValue += t;

            return featureValue;
        }

        private double feature10(List<double> thetas)
        {
            double featureValue = 0.0;

            foreach (double t in thetas)
                featureValue += Math.Abs(t);

            return featureValue;
        }

        private double feature11(List<double> thetas)
        {
            double featureValue = 0.0;

            foreach (double t in thetas)
                featureValue += Math.Pow(t, 2.0);

            return featureValue;
        }

        private double feature12(List<Sketch.Point> points)
        {
            double featureValue = 0.0;

            double deltax = 0.0;
            double deltay = 0.0;
            double deltat = 0.0;

            for (int i = 0; i < points.Count - 1; i++)
            {
                deltax = Math.Pow((double)(points[i + 1].X - points[i].X), 2.0);
                deltay = Math.Pow((double)(points[i + 1].Y - points[i].Y), 2.0);
                deltat = Math.Pow((double)(points[i + 1].Time - points[i].Time), 2.0);

                if ((deltax + deltay) / deltat > featureValue)
                    featureValue = (deltax + deltay) / deltat;
            }

            return featureValue;
        }

        private double feature13(List<Sketch.Point> points)
        {
            return (double)(points[points.Count - 1].Time - points[0].Time);
        }
        #endregion
    }
}
