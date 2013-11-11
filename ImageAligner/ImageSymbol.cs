using System;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using System.Diagnostics;
using Utilities;
using Utilities.Matrix;
using RecognitionTemplates;

namespace ImageAligner
{
    [DebuggerDisplay("{base.Name}")]
    [Serializable]
    public class ImageSymbol : RecognitionTemplate
    {
        #region Member Variables & Constants

        #region PARAMETERS

        private int GRID_SIZE = 48;
        private const double REL_DIST_SCALE_FACTOR = 2.4;
        private const double HAUSDORFF_QUANTILE = 0.94;

        #endregion

        private SymbolInfo _SymbolInfo;

        /// <summary>
        /// Points that have been quantized to the specified grid size
        /// </summary>
        private List<Point> _screenQuantizedPoints;

        /// <summary>
        /// Distance Transform Matrix for fast computation of an 
        /// overlayed symbols distances
        /// </summary>
        private GeneralMatrix _sDTM;

        /// <summary>
        /// Points that have been quantized to the specified grid size
        /// </summary>
        //private List<Point> _polarQuantizedPoints;

        /// <summary>
        /// Distance Transform Matrix for fast computation of an 
        /// overlayed symbols distances
        /// </summary>
        //private GeneralMatrix _pDTM;

        //private int _polarYmax;

        #endregion

        #region Constructors

        public ImageSymbol(Point[] points, Utilities.SymbolInfo info) : base(info)
        {
            if (backtrack.Count == 0)
                FillBacktrack();

            _SymbolInfo = new SymbolInfo(info.User, info.SymbolType, info.SymbolClass);
            Process(points);
        }

        public ImageSymbol(Point[] points, Rectangle bbox, SymbolInfo info)
        {
            if (backtrack.Count == 0)
                FillBacktrack();

            _SymbolInfo = info;
            BoundingPoints corners = new BoundingPoints(bbox);
            Process(points, corners);
        }

        public ImageSymbol(List<Sketch.Substroke> strokes, Rectangle bbox, SymbolInfo info)
        {
            if (backtrack.Count == 0)
                FillBacktrack();

            _SymbolInfo = info;
            BoundingPoints corners = new BoundingPoints(bbox);
            Point[] points = GetPoints(strokes);
            Process(points, corners);
        }

        public ImageSymbol()
        {
            if (backtrack.Count == 0)
                FillBacktrack();
        }

        #endregion

        #region Getters

        private List<Point> QuantizedScreenPoints
        {
            get 
            {
                //List<Point> points = new List<Point>(_screenQuantizedPoints.Length / 2);
                //for (int i = 0; i < _screenQuantizedPoints.Length / 2; i++)
                    //points[i] = new Point(_screenQuantizedPoints[i, 0], _screenQuantizedPoints[i, 1]);

                return _screenQuantizedPoints;
            }
        }

        private GeneralMatrix DistanceTransformMatrixScreen
        {
            get { return _sDTM; }
        }
        
        public override string Name
        {
            get { return _SymbolInfo.SymbolType; }
        }

        #endregion

        #region Processing

        private Point[] GetPoints(List<Sketch.Substroke> strokes)
        {
            List<Point> points = new List<Point>();
            foreach (Sketch.Substroke stroke in strokes)
                foreach (Sketch.Point pt in stroke.Points)
                    points.Add(pt.SysDrawPoint);

            return points.ToArray();
        }

        /// <summary>
        /// Processes the original points to get the "Bitmap" image of the shape
        /// in both screen and polar coordinates
        /// </summary>
        /// <param name="points">Original points in the shape</param>
        private void Process(Point[] points)
        {
            BoundingPoints corners = new BoundingPoints(points);
            Process(points, corners);
        }

        /// <summary>
        /// Processes the original points to get the "Bitmap" image of the shape
        /// in both screen and polar coordinates
        /// </summary>
        /// <param name="points">Original points in the shape</param>
        /// <param name="corners">Boundaries of the shape</param>
        private void Process(Point[] points, BoundingPoints corners)
        {
            _screenQuantizedPoints = QuantizePointsScreen(points, corners);
            
            //List<Point> qPoints = new List<Point>();
            //for (int i = 0; i < _screenQuantizedPoints.Count; i++)
                //qPoints.Add(new Point(_screenQuantizedPoints[i, 1], _screenQuantizedPoints[i, 0]));
            //GeneralMatrix slow = DistanceTransformScreen(_screenQuantizedPoints);

            _sDTM = DistanceTransformMatrix(_screenQuantizedPoints);
            //GeneralMatrix fast = DistanceTransformMatrix(_screenQuantizedPoints);

            //List<PolarPoint> polarPoints = Transform_Screen2Polar(points);
            
            /*GeneralMatrix pointsSlow = new GeneralMatrix(GRID_SIZE, GRID_SIZE, 0.0);
            foreach (Point pt in _screenQuantizedPoints)
                pointsSlow.SetElement(pt.Y, pt.X, 1.0);

            System.IO.StreamWriter writer = new System.IO.StreamWriter("C:/arrays.txt");

            writer.WriteLine("Points Slow");
            writer.WriteLine(pointsSlow.ToString());
            writer.WriteLine(); writer.WriteLine("DTM Slow");
            writer.WriteLine(slow.ToString());
            writer.WriteLine(); writer.WriteLine("DTM Fast");
            writer.WriteLine(_sDTM.ToString());

            writer.Close();

            for (int i = 0; i < slow.RowDimension; i++)
                for (int j = 0; j < slow.ColumnDimension; j++)
                    if (slow.GetElement(i, j) != _sDTM.GetElement(i, j))
                        Console.WriteLine("Diff = " + (slow.GetElement(i, j) - _sDTM.GetElement(i, j)));
            */
        }

        /// <summary>
        /// Determines the coordinates of the "rasterized" points (in screen coordinates)
        /// </summary>
        /// <param name="points">Original screen points</param>
        /// <param name="corners">The corners of the shape's bounding box</param>
        /// <returns>Quantized points for screen coordinate system</returns>
        private List<Point> QuantizePointsScreen(Point[] points, BoundingPoints corners)
        {
            if (points.Length == 0) return new List<Point>();

            double sq_side = Math.Max(corners.Height, corners.Width);
            double step = sq_side / (double)(GRID_SIZE - 1);
            GeneralMatrix mesh = new GeneralMatrix(GRID_SIZE, GRID_SIZE, 0.0);
            float max = (float)step / 2;

            List<PointF> fPoints = new List<PointF>(points.Length);
            fPoints.Add(new PointF(points[0].X, points[0].Y));
            for (int i = 1; i < points.Length; i++)
            {
                //List<PointF> ipt = Interpolate(points[i], points[i - 1], 1, 1);
                List<PointF> ipt = Interpolate(points[i], points[i - 1], max);
                foreach (PointF pt in ipt)
                    fPoints.Add(pt);

                fPoints.Add(points[i]);
            }
            //foreach (Point pt in points)
                //fPoints.Add(pt);

            PointF c = corners.Center;

            // For each points calculate its relative location inside the bounding box
            for (int i = 0; i < fPoints.Count; i++)
            {
                int x_index = (int)Math.Floor(((fPoints[i].X - c.X) / step) + GRID_SIZE / 2);
                int y_index = (int)Math.Floor(((fPoints[i].Y - c.Y) / step) + GRID_SIZE / 2);

                if (x_index < 0)
                    x_index = 0;
                else if (x_index >= GRID_SIZE)
                    x_index = GRID_SIZE - 1;

                if (y_index < 0)
                    y_index = 0;
                else if (y_index >= GRID_SIZE)
                    y_index = GRID_SIZE - 1;

                mesh.SetElement(y_index, x_index, 1.0);
            }

            List<Point> qPoints = new List<Point>();
            // Go throught the entire matrix and create new points 
            // whenever you encounter a value greater than 0
            for (int i = 0; i < GRID_SIZE; i++)
                for (int j = 0; j < GRID_SIZE; j++)
                    if (mesh.GetElement(i, j) > 0.0)
                        qPoints.Add(new Point(j, i));

            return qPoints;
        }

        private List<PointF> Interpolate(Point point, Point point_2, double max)
        {
            double dist = Compute.EuclideanDistance(point, point_2);
            int numPoints = (int)Math.Floor(dist / max);

            List<PointF> points = new List<PointF>(numPoints);

            if (dist > 200.0)
                return points;

            float dx = point_2.X - point.X;
            dx /= (numPoints + 1);
            float dy = point_2.Y - point.Y;
            dy /= (numPoints + 1);

            for (int i = 1; i <= numPoints; i++)
            {
                float x = point.X + dx * i;
                float y = point.Y + dy * i;
                points.Add(new PointF(x, y));
            }

            return points;
        }

        private List<PointF> Interpolate(PointF point, PointF point_2, int numHalves, int currentLayer)
        {
            if (currentLayer >= numHalves)
                return new List<PointF>(new PointF[] { Midpoint(point, point_2) });

            int numPoints = 2;
            for (int i = 1; i <= numHalves - currentLayer; i++)
                numPoints += i * i;

            List<PointF> points = new List<PointF>(numPoints);
            points.Add(point);

            PointF mid = Midpoint(point, point_2);
            
            List<PointF> pts1 = Interpolate(point, mid, numHalves, currentLayer + 1);
            foreach (PointF pt in pts1)
                points.Add(pt);

            points.Add(mid);

            List<PointF> pts2 = Interpolate(mid, point_2, numHalves, currentLayer + 1);
            foreach (PointF pt in pts2)
                points.Add(pt);

            points.Add(point_2);
            
            return points;
        }

        private PointF Midpoint(PointF point1, PointF point2)
        {
            float x = (point2.X + point1.X) / 2;
            float y = (point2.Y + point1.Y) / 2;
            return new PointF(x, y);
        }

        /// <summary>
        /// Determines the coordinates of the "rasterized" points (in polar coordinates)
        /// </summary>
        /// <param name="points">Regular polar points</param>
        /// <returns>Quantized points for polar coordinate system</returns>
        private List<Point> QuantizePointsPolar(List<PolarPoint> points)
        {
            List<Point> qPoints = new List<Point>();
            if (points.Count == 0) return qPoints;

            GeneralMatrix mesh = new GeneralMatrix(GRID_SIZE, GRID_SIZE, 0.0);

            double stepX = 2.0 * Math.PI / GRID_SIZE;
            double stepY = REL_DIST_SCALE_FACTOR / GRID_SIZE;

            // For each points calculate its relative location inside the bounding box
            foreach (PolarPoint pt in points)
            {
                int x_index = (int)Math.Floor(((pt.AngularPosition + Math.PI) / stepX));
                int y_index = (int)Math.Floor((pt.RelativeDistance / stepY));

                if (x_index < 0)
                    x_index = 0;
                else if (x_index >= GRID_SIZE)
                    x_index = GRID_SIZE - 1;

                if (y_index < 0)
                    y_index = 0;
                else if (y_index >= GRID_SIZE)
                    y_index = GRID_SIZE - 1;

                mesh.SetElement(y_index, x_index, 1.0);
            }

            // Go throught the entire matrix and create new points 
            // whenever you encounter a value greater than 0
            for (int i = 0; i < GRID_SIZE; i++)
                for (int j = 0; j < GRID_SIZE; j++)
                    if (mesh.GetElement(i, j) > 0.0)
                        qPoints.Add(new Point(j, i));

            return qPoints;
        }

        /// <summary>
        /// Computes the distance transform matrix for screen coordinates
        /// </summary>
        /// <param name="qPoints">quantized points in screen coordinates</param>
        /// <returns>Distance Transform Matrix</returns>
        private GeneralMatrix DistanceTransformScreen(List<Point> qPoints)
        {
            GeneralMatrix DTM = new GeneralMatrix(GRID_SIZE, GRID_SIZE, double.PositiveInfinity);

            for (int i = 0; i < GRID_SIZE; i++)
            {
                for (int j = 0; j < GRID_SIZE; j++)
                {
                    Point current = new Point(j, i);

                    double mindist = double.PositiveInfinity;

                    foreach (Point pt in qPoints)
                        mindist = Math.Min(mindist, Compute.EuclideanDistance(current, pt));

                    DTM.SetElement(i, j, mindist);
                }
            }

            return DTM;
        }

        private GeneralMatrix DistanceTransformMatrix(List<Point> points)
        {
            GeneralMatrix pMat = new GeneralMatrix(GRID_SIZE, GRID_SIZE, 0.0);
            for (int i = 0; i < points.Count; i++)
                pMat.SetElement(points[i].Y, points[i].X, 1.0);

            return DistanceTransformMatrix(pMat);
        }

        private GeneralMatrix DistanceTransformMatrix(GeneralMatrix pMat)
        {
            GeneralMatrix DTM = new GeneralMatrix(GRID_SIZE, GRID_SIZE, double.PositiveInfinity);

            int start = 0;
            int startRowAbove = 0;

            for (int i = 0; i < GRID_SIZE; i++)
            {
                start = backtrack[startRowAbove];
                for (int j = 0; j < GRID_SIZE; j++)
                {
                    DTM.SetElement(i, j, FindDistanceToNearestPoint(pMat, i, j, ref start));
                    if (j == 0)
                        startRowAbove = start;
                }
            }

            return DTM;
        }

        private double FindDistanceToNearestPoint(GeneralMatrix pMat, int i, int j, ref int start)
        {
            if (pMat.GetElement(i, j) > 0.0)
                return 0.0;

            for (int k = start; k < distances.Length / 2; k++)
                if (PointFoundAtDistance(pMat, i, j, distances[k, 0], distances[k, 1]))
                {
                    start = backtrack[k];
                    return Math.Sqrt((double)(distances[k, 0] * distances[k, 0] + distances[k, 1] * distances[k, 1]));
                }

            return double.MaxValue;
        }

        private bool PointFoundAtDistance(GeneralMatrix points, int i, int j, int d1, int d2)
        {
            int x = i + d1;
            int y = j + d2;
            if (x < points.ColumnDimension && y < points.RowDimension && points.GetElement(x, y) > 0.0)
                return true;

            y = j - d2;
            if (x < points.ColumnDimension && y >= 0 && points.GetElement(x, y) > 0.0)
                return true;

            x = i - d1;
            if (x >= 0 && y >= 0 && points.GetElement(x, y) > 0.0)
                return true;

            y = j + d2;
            if (x >= 0 && y < points.RowDimension && points.GetElement(x, y) > 0.0)
                return true;

            x = i + d2;
            y = j + d1;
            if (x < points.ColumnDimension && y < points.RowDimension && points.GetElement(x, y) > 0.0)
                return true;

            y = j - d1;
            if (x < points.ColumnDimension && y >= 0 && points.GetElement(x, y) > 0.0)
                return true;

            x = i - d2;
            if (x >= 0 && y >= 0 && points.GetElement(x, y) > 0.0)
                return true;

            y = j + d1;
            if (x >= 0 && y < points.RowDimension && points.GetElement(x, y) > 0.0)
                return true;

            return false;
        }

        /// <summary>
        /// Computes the distance transform matrix for polar coordinates.
        /// In polar coordinates this means we have to check for a shift of 
        /// 2 pi
        /// </summary>
        /// <param name="qPoints">quantized points in polar coordinates</param>
        /// <returns>Distance Transform Matrix</returns>
        private GeneralMatrix DistanceTransformPolar(List<Point> qPoints)
        {
            GeneralMatrix DTM = new GeneralMatrix(GRID_SIZE, GRID_SIZE, double.PositiveInfinity);

            for (int i = 0; i < GRID_SIZE; i++)
            {
                for (int j = 0; j < GRID_SIZE; j++)
                {
                    Point current = new Point(j, i);

                    double mindist = double.PositiveInfinity;

                    foreach (Point pt in qPoints)
                    {
                        double distance1 = Compute.EuclideanDistance(current, pt);
                        double dx = (double)GRID_SIZE - Math.Abs(current.X - pt.X);
                        double dy = Math.Abs(current.Y - pt.Y);
                        double distance2 = Math.Sqrt(Math.Pow(dx, 2.0) + Math.Pow(dy, 2.0));
                        double distance = Math.Min(distance1, distance2);
                        mindist = Math.Min(mindist, distance);
                    }

                    DTM.SetElement(i, j, mindist);
                }
            }

            return DTM;
        }

        /// <summary>
        /// Creates polar points from the screen coordinates
        /// </summary>
        /// <param name="points">points in screen coordinates</param>
        /// <returns>Polar points - unquantized</returns>
        private List<PolarPoint> Transform_Screen2Polar(Point[] points)
        {
            if (points.Length == 0)
                return new List<PolarPoint>();

            int sumX = 0;
            int sumY = 0;
            foreach (Point pt in points)
            {
                sumX += pt.X;
                sumY += pt.Y;
            }
            Point center = new Point(sumX / points.Length, sumY / points.Length);

            double avgDistance = 0.0;

            Dictionary<Point, double> distances = new Dictionary<Point, double>(points.Length);
            foreach (Point pt in points)
            {
                double dist = Compute.EuclideanDistance(pt, center);
                if (!distances.ContainsKey(pt))
                    distances.Add(pt, dist);
                avgDistance += dist;
            }

            avgDistance /= points.Length;
            if (avgDistance == 0.0)
            {
                avgDistance += 0.00001;
            }

            List<PolarPoint> polarPoints = new List<PolarPoint>(points.Length);
            foreach (Point pt in points)
            {
                double angPosition = Math.Atan2(pt.Y - center.Y, pt.X - center.X);
                double relDistance;
                if (distances.ContainsKey(pt))
                    relDistance = distances[pt] / avgDistance;
                else
                    relDistance = Compute.EuclideanDistance(pt, center) / avgDistance;
                polarPoints.Add(new PolarPoint(angPosition, relDistance));
            }

            return polarPoints;
        }

        /// <summary>
        /// Find the maximum "Y" value of the polar quantized points
        /// </summary>
        /// <param name="points">polar quantized points</param>
        /// <returns>Max value - used for scaling</returns>
        private int FindYmax(List<Point> points)
        {
            int max = int.MinValue;
            foreach (Point pt in points)
                max = Math.Max(pt.Y, max);

            return max;
        }

        #endregion

        #region Recognition

        public override double Recognize(RecognitionTemplate unknown)
        {
            ImageSymbol u = (ImageSymbol)unknown;
            double modHaus = ModifiedHausdorffDistance(u);
            double haus = HausdorffDistance(u);
            //ImageSymbolResult result = Recognize((ImageSymbol)unknown);
            double score = 1.0 / (1.0 + haus * ImageSymbolResult.HAUS_SCALE + modHaus * ImageSymbolResult.MODHAUS_SCALE);
            return score;
        }

        /// <summary>
        /// Assume no rotation!
        /// </summary>
        /// <param name="templates"></param>
        /// <returns></returns>
        public SortedList<double, ImageSymbolResult> Recognize(List<ImageSymbol> templates)
        {
            SortedList<double, ImageSymbolResult> results = new SortedList<double, ImageSymbolResult>(templates.Count);

            foreach (ImageSymbol symbol in templates)
            {
                ImageSymbolResult r = Recognize(symbol);
                double score = r.Score;
                while (results.ContainsKey(score))
                    score += double.Epsilon;

                results.Add(score, r);
            }

            return results;
        }

        /// <summary>
        /// Calculates each of the image-based metrics between two symbols
        /// </summary>
        /// <param name="template">ImageSymbol to compare 'this' to</param>
        /// <returns>ImageResult for the two symbols</returns>
        public ImageSymbolResult Recognize(ImageSymbol template)
        {
            double haus = HausdorffDistance(template);
            double modHaus = ModifiedHausdorffDistance(template);
            double tanimoto = TanimotoCoefficient(template);
            double yule = YuleCoefficient(template);

            return new ImageSymbolResult(this, template, haus, modHaus, tanimoto, yule);
        }

        /// <summary>
        /// Calculates the Hausdorff Distance (HD) between the two ImageSymbols.
        /// HD is defined as the max pixel distance (in Xth quantile) between the images (max of 2 directed).
        /// </summary>
        /// <param name="template">ImageSymbol to compare 'this' to</param>
        /// <returns>Hausdorff Distance</returns>
        private double HausdorffDistance(ImageSymbol template)
        {
            double hAB = this.DirectedHausdorffDistance(template);
            double hBA = template.DirectedHausdorffDistance(this);

            return Math.Max(hAB, hBA);
        }

        /// <summary>
        /// Calculates the directed Hausdorff Distance between the two ImageSymbols.
        /// Uses a quantile to filter out most extreme outliers
        /// </summary>
        /// <param name="template">ImageSymbol to compare 'this' to</param>
        /// <returns>Directed Hausdorff Distance</returns>
        private double DirectedHausdorffDistance(ImageSymbol template)
        {
            List<double> distances = new List<double>(template._screenQuantizedPoints.Count);
            for (int i = 0; i < template._screenQuantizedPoints.Count; i++)
            {
                int x = template._screenQuantizedPoints[i].X;
                int y = template._screenQuantizedPoints[i].Y;
                if (y < _sDTM.RowDimension && y >= 0 && x >= 0 && x < _sDTM.ColumnDimension)
                {
                    double distan = _sDTM.GetElement(y, x);
                    distances.Add(distan);
                }
            }

            distances.Sort();
            if (distances.Count == 0) return double.PositiveInfinity;

            return distances[(int)Math.Floor(((distances.Count - 1) * HAUSDORFF_QUANTILE))];
        }

        /// <summary>
        /// Calculates the Modified Hausdorff Distance (MHD) between the two ImageSymbols.
        /// MHD is defined as the average pixel distance between the images (max of 2 directed).
        /// </summary>
        /// <param name="template">ImageSymbol to compare 'this' to</param>
        /// <returns>Modified Hausdorff Distance</returns>
        private double ModifiedHausdorffDistance(ImageSymbol template)
        {
            double distan1 = 0.0;
            for (int i = 0; i < _screenQuantizedPoints.Count; i++)
            {
                int x = _screenQuantizedPoints[i].X;
                int y = _screenQuantizedPoints[i].Y;
                if (y >= 0 && x >= 0
                    && y < template.DistanceTransformMatrixScreen.RowDimension
                    && x < template.DistanceTransformMatrixScreen.ColumnDimension)
                    distan1 += template.DistanceTransformMatrixScreen.GetElement(y, x);
            }

            double AB = distan1 / (_screenQuantizedPoints.Count);

            double distan2 = 0.0;
            for (int i = 0; i < template._screenQuantizedPoints.Count; i++)
            {
                int x = template._screenQuantizedPoints[i].X;
                int y = template._screenQuantizedPoints[i].Y;
                if (y >= 0 && x >= 0 && y < _sDTM.RowDimension && x < _sDTM.ColumnDimension)
                    distan2 += _sDTM.GetElement(y, x);
            }

            double BA = distan2 / (template._screenQuantizedPoints.Count);

            return Math.Max(AB, BA);
        }

        /// <summary>
        /// Calculates the Tanimoto coefficient between the two ImageSymbols
        /// </summary>
        /// <param name="template">ImageSymbol to compare 'this' to</param>
        /// <returns>Tanimoto Coefficient</returns>
        private double TanimotoCoefficient(ImageSymbol template)
        {
            int a, b, c, d;
            a = b = c = d = 0;
            double E = 1.0 / 15.0 * Math.Sqrt(Math.Pow(GRID_SIZE, 2.0) + Math.Pow(GRID_SIZE, 2.0));

            for (int i = 0; i < _sDTM.ColumnDimension; i++)
            {
                for (int j = 0; j < _sDTM.RowDimension; j++)
                {
                    if (_sDTM.GetElement(i, j) < E)
                        a++;

                    if (template.DistanceTransformMatrixScreen.GetElement(i, j) < E)
                        b++;

                    if (_sDTM.GetElement(i, j) < E && template.DistanceTransformMatrixScreen.GetElement(i, j) < E)
                        c++;

                    if (_sDTM.GetElement(i, j) >= E && template.DistanceTransformMatrixScreen.GetElement(i, j) >= E)
                        d++;
                }
            }

            double T1 = (double)c / (double)(a + b - c);
            double T0 = (double)d / (double)(a + b - 2 * c + d);
            double p = (double)(a + b) / 2.0 / _sDTM.ColumnDimension / _sDTM.RowDimension;
            double alpha = (3.0 - p) / 4.0;

            return 1.0 - (alpha * T1 + (1.0 - alpha) * T0);
        }

        /// <summary>
        /// Calculates the Yule coefficient between the two ImageSymbols
        /// </summary>
        /// <param name="template">ImageSymbol to compare 'this' to</param>
        /// <returns>Yule Coefficient</returns>
        private double YuleCoefficient(ImageSymbol template)
        {
            int n10, n01, n11, n00;
            n10 = n01 = n11 = n00 = 0;
            double E = 1.0 / 15.0 * Math.Sqrt(Math.Pow(GRID_SIZE, 2.0) + Math.Pow(GRID_SIZE, 2.0));

            for (int i = 0; i < _sDTM.ColumnDimension; i++)
            {
                for (int j = 0; j < _sDTM.RowDimension; j++)
                {
                    if (_sDTM.GetElement(i, j) < E && template.DistanceTransformMatrixScreen.GetElement(i, j) >= E)
                        n10++;

                    if (_sDTM.GetElement(i, j) >= E && template.DistanceTransformMatrixScreen.GetElement(i, j) < E)
                        n01++;

                    if (_sDTM.GetElement(i, j) < E && template.DistanceTransformMatrixScreen.GetElement(i, j) < E)
                        n11++;

                    if (_sDTM.GetElement(i, j) >= E && template.DistanceTransformMatrixScreen.GetElement(i, j) >= E)
                        n00++;
                }
            }

            return 1.0 - (double)(n11 * n00 - n10 * n01) / (double)(n11 * n00 + n10 * n01);
        }

        #endregion

        #region Make Actual Bitmap Image of Symbol

        public Bitmap MakeBitmap()
        {
            Bitmap map = new Bitmap(GRID_SIZE, GRID_SIZE);

            Graphics DrawSurface = Graphics.FromImage(map);

            DrawSurface.Clear(Color.White);
            Pen p = new Pen(Color.Black);
            Brush b = Brushes.Black;

            for (int i = 0; i < _screenQuantizedPoints.Count; i++)
            {
                DrawSurface.FillEllipse(b, _screenQuantizedPoints[i].X, _screenQuantizedPoints[i].Y, 5, 5);
                
                //map.SetPixel(_screenQuantizedPoints[i, 0], _screenQuantizedPoints[i, 1], Color.Black);
            }

            

            return map;
        }

        #endregion

        #region Serialization



        #endregion

        #region STATIC CONSTANTS

        private static Dictionary<int, int> backtrack = new Dictionary<int, int>();

        private void FillBacktrack()
        {
            backtrack.Add(0, 0);
            backtrack.Add(1, 0);
            backtrack.Add(2, 0);
            backtrack.Add(3, 1);
            backtrack.Add(4, 1);
            backtrack.Add(5, 1);
            backtrack.Add(6, 3);
            backtrack.Add(7, 3);
            backtrack.Add(8, 3);
            backtrack.Add(9, 6);
            backtrack.Add(10, 6);
            backtrack.Add(11, 6);
            backtrack.Add(12, 6);
            backtrack.Add(13, 9);
            backtrack.Add(14, 9);
            backtrack.Add(15, 9);
            backtrack.Add(16, 9);
            backtrack.Add(17, 9);
            backtrack.Add(18, 9);
            backtrack.Add(19, 13);
            backtrack.Add(20, 13);
            backtrack.Add(21, 13);
            backtrack.Add(22, 13);
            backtrack.Add(23, 13);
            backtrack.Add(24, 19);
            backtrack.Add(25, 19);
            backtrack.Add(26, 19);
            backtrack.Add(27, 19);
            backtrack.Add(28, 19);
            backtrack.Add(29, 19);
            backtrack.Add(30, 19);
            backtrack.Add(31, 24);
            backtrack.Add(32, 24);
            backtrack.Add(33, 24);
            backtrack.Add(34, 24);
            backtrack.Add(35, 24);
            backtrack.Add(36, 24);
            backtrack.Add(37, 24);
            backtrack.Add(38, 24);
            backtrack.Add(39, 31);
            backtrack.Add(40, 31);
            backtrack.Add(41, 31);
            backtrack.Add(42, 31);
            backtrack.Add(43, 31);
            backtrack.Add(44, 31);
            backtrack.Add(45, 31);
            backtrack.Add(46, 31);
            backtrack.Add(47, 39);
            backtrack.Add(48, 39);
            backtrack.Add(49, 39);
            backtrack.Add(50, 39);
            backtrack.Add(51, 39);
            backtrack.Add(52, 39);
            backtrack.Add(53, 39);
            backtrack.Add(54, 39);
            backtrack.Add(55, 39);
            backtrack.Add(56, 47);
            backtrack.Add(57, 47);
            backtrack.Add(58, 47);
            backtrack.Add(59, 47);
            backtrack.Add(60, 47);
            backtrack.Add(61, 47);
            backtrack.Add(62, 47);
            backtrack.Add(63, 47);
            backtrack.Add(64, 47);
            backtrack.Add(65, 56);
            backtrack.Add(66, 56);
            backtrack.Add(67, 56);
            backtrack.Add(68, 56);
            backtrack.Add(69, 56);
            backtrack.Add(70, 56);
            backtrack.Add(71, 56);
            backtrack.Add(72, 56);
            backtrack.Add(73, 56);
            backtrack.Add(74, 56);
            backtrack.Add(75, 56);
            backtrack.Add(76, 65);
            backtrack.Add(77, 65);
            backtrack.Add(78, 65);
            backtrack.Add(79, 65);
            backtrack.Add(80, 65);
            backtrack.Add(81, 65);
            backtrack.Add(82, 65);
            backtrack.Add(83, 65);
            backtrack.Add(84, 65);
            backtrack.Add(85, 65);
            backtrack.Add(86, 65);
            backtrack.Add(87, 65);
            backtrack.Add(88, 76);
            backtrack.Add(89, 76);
            backtrack.Add(90, 76);
            backtrack.Add(91, 76);
            backtrack.Add(92, 76);
            backtrack.Add(93, 76);
            backtrack.Add(94, 76);
            backtrack.Add(95, 76);
            backtrack.Add(96, 76);
            backtrack.Add(97, 76);
            backtrack.Add(98, 76);
            backtrack.Add(99, 76);
            backtrack.Add(100, 88);
            backtrack.Add(101, 88);
            backtrack.Add(102, 88);
            backtrack.Add(103, 88);
            backtrack.Add(104, 88);
            backtrack.Add(105, 88);
            backtrack.Add(106, 88);
            backtrack.Add(107, 88);
            backtrack.Add(108, 88);
            backtrack.Add(109, 88);
            backtrack.Add(110, 88);
            backtrack.Add(111, 88);
            backtrack.Add(112, 88);
            backtrack.Add(113, 100);
            backtrack.Add(114, 100);
            backtrack.Add(115, 100);
            backtrack.Add(116, 100);
            backtrack.Add(117, 100);
            backtrack.Add(118, 100);
            backtrack.Add(119, 100);
            backtrack.Add(120, 100);
            backtrack.Add(121, 100);
            backtrack.Add(122, 100);
            backtrack.Add(123, 100);
            backtrack.Add(124, 100);
            backtrack.Add(125, 100);
            backtrack.Add(126, 113);
            backtrack.Add(127, 113);
            backtrack.Add(128, 113);
            backtrack.Add(129, 113);
            backtrack.Add(130, 113);
            backtrack.Add(131, 113);
            backtrack.Add(132, 113);
            backtrack.Add(133, 113);
            backtrack.Add(134, 113);
            backtrack.Add(135, 113);
            backtrack.Add(136, 113);
            backtrack.Add(137, 113);
            backtrack.Add(138, 113);
            backtrack.Add(139, 113);
            backtrack.Add(140, 113);
            backtrack.Add(141, 126);
            backtrack.Add(142, 126);
            backtrack.Add(143, 126);
            backtrack.Add(144, 126);
            backtrack.Add(145, 126);
            backtrack.Add(146, 126);
            backtrack.Add(147, 126);
            backtrack.Add(148, 126);
            backtrack.Add(149, 126);
            backtrack.Add(150, 126);
            backtrack.Add(151, 126);
            backtrack.Add(152, 126);
            backtrack.Add(153, 126);
            backtrack.Add(154, 126);
            backtrack.Add(155, 126);
            backtrack.Add(156, 126);
            backtrack.Add(157, 141);
            backtrack.Add(158, 141);
            backtrack.Add(159, 141);
            backtrack.Add(160, 141);
            backtrack.Add(161, 141);
            backtrack.Add(162, 141);
            backtrack.Add(163, 141);
            backtrack.Add(164, 141);
            backtrack.Add(165, 141);
            backtrack.Add(166, 141);
            backtrack.Add(167, 141);
            backtrack.Add(168, 141);
            backtrack.Add(169, 141);
            backtrack.Add(170, 141);
            backtrack.Add(171, 141);
            backtrack.Add(172, 141);
            backtrack.Add(173, 157);
            backtrack.Add(174, 157);
            backtrack.Add(175, 157);
            backtrack.Add(176, 157);
            backtrack.Add(177, 157);
            backtrack.Add(178, 157);
            backtrack.Add(179, 157);
            backtrack.Add(180, 157);
            backtrack.Add(181, 157);
            backtrack.Add(182, 157);
            backtrack.Add(183, 157);
            backtrack.Add(184, 157);
            backtrack.Add(185, 157);
            backtrack.Add(186, 157);
            backtrack.Add(187, 157);
            backtrack.Add(188, 157);
            backtrack.Add(189, 173);
            backtrack.Add(190, 173);
            backtrack.Add(191, 173);
            backtrack.Add(192, 173);
            backtrack.Add(193, 173);
            backtrack.Add(194, 173);
            backtrack.Add(195, 173);
            backtrack.Add(196, 173);
            backtrack.Add(197, 173);
            backtrack.Add(198, 173);
            backtrack.Add(199, 173);
            backtrack.Add(200, 173);
            backtrack.Add(201, 173);
            backtrack.Add(202, 173);
            backtrack.Add(203, 173);
            backtrack.Add(204, 173);
            backtrack.Add(205, 173);
            backtrack.Add(206, 173);
            backtrack.Add(207, 173);
            backtrack.Add(208, 189);
            backtrack.Add(209, 189);
            backtrack.Add(210, 189);
            backtrack.Add(211, 189);
            backtrack.Add(212, 189);
            backtrack.Add(213, 189);
            backtrack.Add(214, 189);
            backtrack.Add(215, 189);
            backtrack.Add(216, 189);
            backtrack.Add(217, 189);
            backtrack.Add(218, 189);
            backtrack.Add(219, 189);
            backtrack.Add(220, 189);
            backtrack.Add(221, 189);
            backtrack.Add(222, 189);
            backtrack.Add(223, 189);
            backtrack.Add(224, 189);
            backtrack.Add(225, 189);
            backtrack.Add(226, 208);
            backtrack.Add(227, 208);
            backtrack.Add(228, 208);
            backtrack.Add(229, 208);
            backtrack.Add(230, 208);
            backtrack.Add(231, 208);
            backtrack.Add(232, 208);
            backtrack.Add(233, 208);
            backtrack.Add(234, 208);
            backtrack.Add(235, 208);
            backtrack.Add(236, 208);
            backtrack.Add(237, 208);
            backtrack.Add(238, 208);
            backtrack.Add(239, 208);
            backtrack.Add(240, 208);
            backtrack.Add(241, 208);
            backtrack.Add(242, 208);
            backtrack.Add(243, 208);
            backtrack.Add(244, 226);
            backtrack.Add(245, 226);
            backtrack.Add(246, 226);
            backtrack.Add(247, 226);
            backtrack.Add(248, 226);
            backtrack.Add(249, 226);
            backtrack.Add(250, 226);
            backtrack.Add(251, 226);
            backtrack.Add(252, 226);
            backtrack.Add(253, 226);
            backtrack.Add(254, 226);
            backtrack.Add(255, 226);
            backtrack.Add(256, 226);
            backtrack.Add(257, 226);
            backtrack.Add(258, 226);
            backtrack.Add(259, 226);
            backtrack.Add(260, 226);
            backtrack.Add(261, 226);
            backtrack.Add(262, 226);
            backtrack.Add(263, 226);
            backtrack.Add(264, 244);
            backtrack.Add(265, 244);
            backtrack.Add(266, 244);
            backtrack.Add(267, 244);
            backtrack.Add(268, 244);
            backtrack.Add(269, 244);
            backtrack.Add(270, 244);
            backtrack.Add(271, 244);
            backtrack.Add(272, 244);
            backtrack.Add(273, 244);
            backtrack.Add(274, 244);
            backtrack.Add(275, 244);
            backtrack.Add(276, 244);
            backtrack.Add(277, 244);
            backtrack.Add(278, 244);
            backtrack.Add(279, 244);
            backtrack.Add(280, 244);
            backtrack.Add(281, 244);
            backtrack.Add(282, 244);
            backtrack.Add(283, 244);
            backtrack.Add(284, 244);
            backtrack.Add(285, 244);
            backtrack.Add(286, 264);
            backtrack.Add(287, 264);
            backtrack.Add(288, 264);
            backtrack.Add(289, 264);
            backtrack.Add(290, 264);
            backtrack.Add(291, 264);
            backtrack.Add(292, 264);
            backtrack.Add(293, 264);
            backtrack.Add(294, 264);
            backtrack.Add(295, 264);
            backtrack.Add(296, 264);
            backtrack.Add(297, 264);
            backtrack.Add(298, 264);
            backtrack.Add(299, 264);
            backtrack.Add(300, 264);
            backtrack.Add(301, 264);
            backtrack.Add(302, 264);
            backtrack.Add(303, 264);
            backtrack.Add(304, 264);
            backtrack.Add(305, 264);
            backtrack.Add(306, 264);
            backtrack.Add(307, 264);
            backtrack.Add(308, 264);
            backtrack.Add(309, 286);
            backtrack.Add(310, 286);
            backtrack.Add(311, 286);
            backtrack.Add(312, 286);
            backtrack.Add(313, 286);
            backtrack.Add(314, 286);
            backtrack.Add(315, 286);
            backtrack.Add(316, 286);
            backtrack.Add(317, 286);
            backtrack.Add(318, 286);
            backtrack.Add(319, 286);
            backtrack.Add(320, 286);
            backtrack.Add(321, 286);
            backtrack.Add(322, 286);
            backtrack.Add(323, 286);
            backtrack.Add(324, 286);
            backtrack.Add(325, 286);
            backtrack.Add(326, 286);
            backtrack.Add(327, 286);
            backtrack.Add(328, 286);
            backtrack.Add(329, 286);
            backtrack.Add(330, 309);
            backtrack.Add(331, 309);
            backtrack.Add(332, 309);
            backtrack.Add(333, 309);
            backtrack.Add(334, 309);
            backtrack.Add(335, 309);
            backtrack.Add(336, 309);
            backtrack.Add(337, 309);
            backtrack.Add(338, 309);
            backtrack.Add(339, 309);
            backtrack.Add(340, 309);
            backtrack.Add(341, 309);
            backtrack.Add(342, 309);
            backtrack.Add(343, 309);
            backtrack.Add(344, 309);
            backtrack.Add(345, 309);
            backtrack.Add(346, 309);
            backtrack.Add(347, 309);
            backtrack.Add(348, 309);
            backtrack.Add(349, 309);
            backtrack.Add(350, 309);
            backtrack.Add(351, 309);
            backtrack.Add(352, 330);
            backtrack.Add(353, 330);
            backtrack.Add(354, 330);
            backtrack.Add(355, 330);
            backtrack.Add(356, 330);
            backtrack.Add(357, 330);
            backtrack.Add(358, 330);
            backtrack.Add(359, 330);
            backtrack.Add(360, 330);
            backtrack.Add(361, 330);
            backtrack.Add(362, 330);
            backtrack.Add(363, 330);
            backtrack.Add(364, 330);
            backtrack.Add(365, 330);
            backtrack.Add(366, 330);
            backtrack.Add(367, 330);
            backtrack.Add(368, 330);
            backtrack.Add(369, 330);
            backtrack.Add(370, 330);
            backtrack.Add(371, 330);
            backtrack.Add(372, 330);
            backtrack.Add(373, 330);
            backtrack.Add(374, 330);
            backtrack.Add(375, 330);
            backtrack.Add(376, 330);
            backtrack.Add(377, 352);
            backtrack.Add(378, 352);
            backtrack.Add(379, 352);
            backtrack.Add(380, 352);
            backtrack.Add(381, 352);
            backtrack.Add(382, 352);
            backtrack.Add(383, 352);
            backtrack.Add(384, 352);
            backtrack.Add(385, 352);
            backtrack.Add(386, 352);
            backtrack.Add(387, 352);
            backtrack.Add(388, 352);
            backtrack.Add(389, 352);
            backtrack.Add(390, 352);
            backtrack.Add(391, 352);
            backtrack.Add(392, 352);
            backtrack.Add(393, 352);
            backtrack.Add(394, 352);
            backtrack.Add(395, 352);
            backtrack.Add(396, 352);
            backtrack.Add(397, 352);
            backtrack.Add(398, 352);
            backtrack.Add(399, 352);
            backtrack.Add(400, 352);
            backtrack.Add(401, 377);
            backtrack.Add(402, 377);
            backtrack.Add(403, 377);
            backtrack.Add(404, 377);
            backtrack.Add(405, 377);
            backtrack.Add(406, 377);
            backtrack.Add(407, 377);
            backtrack.Add(408, 377);
            backtrack.Add(409, 377);
            backtrack.Add(410, 377);
            backtrack.Add(411, 377);
            backtrack.Add(412, 377);
            backtrack.Add(413, 377);
            backtrack.Add(414, 377);
            backtrack.Add(415, 377);
            backtrack.Add(416, 377);
            backtrack.Add(417, 377);
            backtrack.Add(418, 377);
            backtrack.Add(419, 377);
            backtrack.Add(420, 377);
            backtrack.Add(421, 377);
            backtrack.Add(422, 377);
            backtrack.Add(423, 377);
            backtrack.Add(424, 377);
            backtrack.Add(425, 377);
            backtrack.Add(426, 377);
            backtrack.Add(427, 377);
            backtrack.Add(428, 401);
            backtrack.Add(429, 401);
            backtrack.Add(430, 401);
            backtrack.Add(431, 401);
            backtrack.Add(432, 401);
            backtrack.Add(433, 401);
            backtrack.Add(434, 401);
            backtrack.Add(435, 401);
            backtrack.Add(436, 401);
            backtrack.Add(437, 401);
            backtrack.Add(438, 401);
            backtrack.Add(439, 401);
            backtrack.Add(440, 401);
            backtrack.Add(441, 401);
            backtrack.Add(442, 401);
            backtrack.Add(443, 401);
            backtrack.Add(444, 401);
            backtrack.Add(445, 401);
            backtrack.Add(446, 401);
            backtrack.Add(447, 401);
            backtrack.Add(448, 401);
            backtrack.Add(449, 401);
            backtrack.Add(450, 401);
            backtrack.Add(451, 401);
            backtrack.Add(452, 401);
            backtrack.Add(453, 401);
            backtrack.Add(454, 428);
            backtrack.Add(455, 428);
            backtrack.Add(456, 428);
            backtrack.Add(457, 428);
            backtrack.Add(458, 428);
            backtrack.Add(459, 428);
            backtrack.Add(460, 428);
            backtrack.Add(461, 428);
            backtrack.Add(462, 428);
            backtrack.Add(463, 428);
            backtrack.Add(464, 428);
            backtrack.Add(465, 428);
            backtrack.Add(466, 428);
            backtrack.Add(467, 428);
            backtrack.Add(468, 428);
            backtrack.Add(469, 428);
            backtrack.Add(470, 428);
            backtrack.Add(471, 428);
            backtrack.Add(472, 428);
            backtrack.Add(473, 428);
            backtrack.Add(474, 428);
            backtrack.Add(475, 428);
            backtrack.Add(476, 428);
            backtrack.Add(477, 428);
            backtrack.Add(478, 428);
        }
        

        private static int[,] distances = new int[479, 2] {
                        { 0, 0 },
                        { 1, 0 },
                        { 1, 1 },
                        { 2, 0 },
                        { 2, 1 },
                        { 2, 2 },
                        { 3, 0 },
                        { 3, 1 },
                        { 3, 2 },
                        { 4, 0 },
                        { 4, 1 },
                        { 3, 3 },
                        { 4, 2 },
                        { 5, 0 },
                        { 4, 3 },
                        { 5, 1 },
                        { 5, 2 },
                        { 4, 4 },
                        { 5, 3 },
                        { 6, 0 },
                        { 6, 1 },
                        { 6, 2 },
                        { 5, 4 },
                        { 6, 3 },
                        { 7, 0 },
                        { 7, 1 },
                        { 5, 5 },
                        { 6, 4 },
                        { 7, 2 },
                        { 7, 3 },
                        { 6, 5 },
                        { 8, 0 },
                        { 8, 1 },
                        { 7, 4 },
                        { 8, 2 },
                        { 6, 6 },
                        { 8, 3 },
                        { 7, 5 },
                        { 8, 4 },
                        { 9, 0 },
                        { 9, 1 },
                        { 9, 2 },
                        { 7, 6 },
                        { 8, 5 },
                        { 9, 3 },
                        { 9, 4 },
                        { 7, 7 },
                        { 10, 0 },
                        { 8, 6 },
                        { 10, 1 },
                        { 10, 2 },
                        { 9, 5 },
                        { 10, 3 },
                        { 8, 7 },
                        { 10, 4 },
                        { 9, 6 },
                        { 11, 0 },
                        { 11, 1 },
                        { 11, 2 },
                        { 10, 5 },
                        { 8, 8 },
                        { 11, 3 },
                        { 9, 7 },
                        { 10, 6 },
                        { 11, 4 },
                        { 12, 0 },
                        { 12, 1 },
                        { 9, 8 },
                        { 11, 5 },
                        { 12, 2 },
                        { 10, 7 },
                        { 12, 3 },
                        { 11, 6 },
                        { 12, 4 },
                        { 9, 9 },
                        { 10, 8 },
                        { 13, 0 },
                        { 12, 5 },
                        { 13, 1 },
                        { 11, 7 },
                        { 13, 2 },
                        { 13, 3 },
                        { 12, 6 },
                        { 10, 9 },
                        { 13, 4 },
                        { 11, 8 },
                        { 12, 7 },
                        { 13, 5 },
                        { 14, 0 },
                        { 14, 1 },
                        { 14, 2 },
                        { 10, 10 },
                        { 11, 9 },
                        { 14, 3 },
                        { 13, 6 },
                        { 12, 8 },
                        { 14, 4 },
                        { 13, 7 },
                        { 14, 5 },
                        { 11, 10 },
                        { 15, 0 },
                        { 12, 9 },
                        { 15, 1 },
                        { 15, 2 },
                        { 14, 6 },
                        { 13, 8 },
                        { 15, 3 },
                        { 15, 4 },
                        { 11, 11 },
                        { 12, 10 },
                        { 14, 7 },
                        { 15, 5 },
                        { 13, 9 },
                        { 16, 0 },
                        { 16, 1 },
                        { 16, 2 },
                        { 14, 8 },
                        { 15, 6 },
                        { 16, 3 },
                        { 12, 11 },
                        { 13, 10 },
                        { 16, 4 },
                        { 15, 7 },
                        { 14, 9 },
                        { 16, 5 },
                        { 12, 12 },
                        { 17, 0 },
                        { 15, 8 },
                        { 17, 1 },
                        { 13, 11 },
                        { 16, 6 },
                        { 17, 2 },
                        { 14, 10 },
                        { 17, 3 },
                        { 17, 4 },
                        { 16, 7 },
                        { 15, 9 },
                        { 13, 12 },
                        { 17, 5 },
                        { 14, 11 },
                        { 16, 8 },
                        { 18, 0 },
                        { 18, 1 },
                        { 17, 6 },
                        { 15, 10 },
                        { 18, 2 },
                        { 18, 3 },
                        { 16, 9 },
                        { 17, 7 },
                        { 13, 13 },
                        { 18, 4 },
                        { 14, 12 },
                        { 15, 11 },
                        { 18, 5 },
                        { 17, 8 },
                        { 16, 10 },
                        { 18, 6 },
                        { 19, 0 },
                        { 19, 1 },
                        { 19, 2 },
                        { 14, 13 },
                        { 15, 12 },
                        { 19, 3 },
                        { 17, 9 },
                        { 18, 7 },
                        { 19, 4 },
                        { 16, 11 },
                        { 19, 5 },
                        { 18, 8 },
                        { 17, 10 },
                        { 14, 14 },
                        { 15, 13 },
                        { 19, 6 },
                        { 20, 0 },
                        { 16, 12 },
                        { 20, 1 },
                        { 20, 2 },
                        { 18, 9 },
                        { 20, 3 },
                        { 19, 7 },
                        { 17, 11 },
                        { 20, 4 },
                        { 15, 14 },
                        { 18, 10 },
                        { 20, 5 },
                        { 19, 8 },
                        { 16, 13 },
                        { 17, 12 },
                        { 20, 6 },
                        { 21, 0 },
                        { 21, 1 },
                        { 19, 9 },
                        { 21, 2 },
                        { 18, 11 },
                        { 20, 7 },
                        { 21, 3 },
                        { 15, 15 },
                        { 16, 14 },
                        { 21, 4 },
                        { 17, 13 },
                        { 19, 10 },
                        { 20, 8 },
                        { 21, 5 },
                        { 18, 12 },
                        { 21, 6 },
                        { 20, 9 },
                        { 16, 15 },
                        { 19, 11 },
                        { 22, 0 },
                        { 22, 1 },
                        { 17, 14 },
                        { 22, 2 },
                        { 21, 7 },
                        { 22, 3 },
                        { 18, 13 },
                        { 22, 4 },
                        { 20, 10 },
                        { 21, 8 },
                        { 19, 12 },
                        { 22, 5 },
                        { 16, 16 },
                        { 17, 15 },
                        { 22, 6 },
                        { 18, 14 },
                        { 20, 11 },
                        { 21, 9 },
                        { 23, 0 },
                        { 23, 1 },
                        { 19, 13 },
                        { 23, 2 },
                        { 22, 7 },
                        { 23, 3 },
                        { 21, 10 },
                        { 20, 12 },
                        { 23, 4 },
                        { 17, 16 },
                        { 22, 8 },
                        { 18, 15 },
                        { 23, 5 },
                        { 19, 14 },
                        { 21, 11 },
                        { 23, 6 },
                        { 22, 9 },
                        { 20, 13 },
                        { 24, 0 },
                        { 24, 1 },
                        { 23, 7 },
                        { 17, 17 },
                        { 24, 2 },
                        { 18, 16 },
                        { 22, 10 },
                        { 24, 3 },
                        { 21, 12 },
                        { 19, 15 },
                        { 24, 4 },
                        { 23, 8 },
                        { 20, 14 },
                        { 24, 5 },
                        { 22, 11 },
                        { 23, 9 },
                        { 21, 13 },
                        { 24, 6 },
                        { 18, 17 },
                        { 19, 16 },
                        { 25, 0 },
                        { 24, 7 },
                        { 20, 15 },
                        { 25, 1 },
                        { 22, 12 },
                        { 25, 2 },
                        { 23, 10 },
                        { 25, 3 },
                        { 21, 14 },
                        { 24, 8 },
                        { 25, 4 },
                        { 18, 18 },
                        { 25, 5 },
                        { 23, 11 },
                        { 19, 17 },
                        { 22, 13 },
                        { 20, 16 },
                        { 24, 9 },
                        { 25, 6 },
                        { 21, 15 },
                        { 23, 12 },
                        { 25, 7 },
                        { 26, 0 },
                        { 24, 10 },
                        { 26, 1 },
                        { 26, 2 },
                        { 22, 14 },
                        { 26, 3 },
                        { 19, 18 },
                        { 25, 8 },
                        { 20, 17 },
                        { 26, 4 },
                        { 24, 11 },
                        { 21, 16 },
                        { 23, 13 },
                        { 26, 5 },
                        { 25, 9 },
                        { 22, 15 },
                        { 26, 6 },
                        { 24, 12 },
                        { 19, 19 },
                        { 20, 18 },
                        { 26, 7 },
                        { 25, 10 },
                        { 23, 14 },
                        { 27, 0 },
                        { 27, 1 },
                        { 21, 17 },
                        { 27, 2 },
                        { 27, 3 },
                        { 26, 8 },
                        { 22, 16 },
                        { 27, 4 },
                        { 24, 13 },
                        { 25, 11 },
                        { 27, 5 },
                        { 23, 15 },
                        { 26, 9 },
                        { 20, 19 },
                        { 27, 6 },
                        { 21, 18 },
                        { 25, 12 },
                        { 24, 14 },
                        { 22, 17 },
                        { 26, 10 },
                        { 27, 7 },
                        { 28, 0 },
                        { 28, 1 },
                        { 23, 16 },
                        { 28, 2 },
                        { 28, 3 },
                        { 27, 8 },
                        { 25, 13 },
                        { 26, 11 },
                        { 28, 4 },
                        { 20, 20 },
                        { 24, 15 },
                        { 21, 19 },
                        { 22, 18 },
                        { 28, 5 },
                        { 27, 9 },
                        { 23, 17 },
                        { 28, 6 },
                        { 26, 12 },
                        { 25, 14 },
                        { 27, 10 },
                        { 24, 16 },
                        { 28, 7 },
                        { 29, 0 },
                        { 21, 20 },
                        { 29, 1 },
                        { 29, 2 },
                        { 26, 13 },
                        { 22, 19 },
                        { 28, 8 },
                        { 29, 3 },
                        { 27, 11 },
                        { 25, 15 },
                        { 23, 18 },
                        { 29, 4 },
                        { 28, 9 },
                        { 24, 17 },
                        { 29, 5 },
                        { 26, 14 },
                        { 27, 12 },
                        { 29, 6 },
                        { 25, 16 },
                        { 21, 21 },
                        { 28, 10 },
                        { 22, 20 },
                        { 29, 7 },
                        { 23, 19 },
                        { 27, 13 },
                        { 30, 0 },
                        { 24, 18 },
                        { 30, 1 },
                        { 26, 15 },
                        { 30, 2 },
                        { 29, 8 },
                        { 28, 11 },
                        { 30, 3 },
                        { 25, 17 },
                        { 30, 4 },
                        { 29, 9 },
                        { 30, 5 },
                        { 27, 14 },
                        { 22, 21 },
                        { 28, 12 },
                        { 23, 20 },
                        { 26, 16 },
                        { 30, 6 },
                        { 24, 19 },
                        { 29, 10 },
                        { 30, 7 },
                        { 25, 18 },
                        { 28, 13 },
                        { 27, 15 },
                        { 31, 0 },
                        { 31, 1 },
                        { 29, 11 },
                        { 30, 8 },
                        { 31, 2 },
                        { 26, 17 },
                        { 22, 22 },
                        { 31, 3 },
                        { 23, 21 },
                        { 24, 20 },
                        { 31, 4 },
                        { 28, 14 },
                        { 30, 9 },
                        { 29, 12 },
                        { 27, 16 },
                        { 31, 5 },
                        { 25, 19 },
                        { 31, 6 },
                        { 30, 10 },
                        { 26, 18 },
                        { 28, 15 },
                        { 31, 7 },
                        { 29, 13 },
                        { 23, 22 },
                        { 24, 21 },
                        { 27, 17 },
                        { 30, 11 },
                        { 32, 0 },
                        { 32, 1 },
                        { 31, 8 },
                        { 25, 20 },
                        { 32, 2 },
                        { 32, 3 },
                        { 29, 14 },
                        { 26, 19 },
                        { 32, 4 },
                        { 28, 16 },
                        { 31, 9 },
                        { 30, 12 },
                        { 32, 5 },
                        { 27, 18 },
                        { 23, 23 },
                        { 32, 6 },
                        { 24, 22 },
                        { 31, 10 },
                        { 29, 15 },
                        { 25, 21 },
                        { 30, 13 },
                        { 32, 7 },
                        { 28, 17 },
                        { 26, 20 },
                        { 31, 11 },
                        { 32, 8 },
                        { 33, 0 },
                        { 33, 1 },
                        { 27, 19 },
                        { 33, 2 },
                        { 30, 14 },
                        { 29, 16 },
                        { 33, 3 },
                        { 33, 4 },
                        { 32, 9 },
                        { 31, 12 },
                        { 24, 23 },
                        { 28, 18 },
                        { 25, 22 },
                        { 33, 5 },
                        { 26, 21 },
                        { 32, 10 },
                        { 33, 6 },
                        { 30, 15 },
                        { 27, 20 },
                        { 31, 13 },
                        { 29, 17 },
                        { 33, 7 },
                        { 32, 11 },
                        { 28, 19 },
                        { 24, 24 }  };

        #endregion
    }

    [Serializable]
    class BoundingPoints
    {
        #region Member Variables

        Point _topLeft;
        Point _topRight;
        Point _bottomLeft;
        Point _bottomRight;

        #endregion

        #region Constructors

        public BoundingPoints(Point[] points)
        {
            FindBoundaries(points);
        }

        public BoundingPoints(Rectangle bbox)
        {
            FindBoundaries(bbox);
        }

        #endregion

        private void FindBoundaries(Point[] points)
        {
            Rectangle bbox = Compute.BoundingBox(points);
            FindBoundaries(bbox);
        }

        private void FindBoundaries(Rectangle bbox)
        {
            _topLeft = new Point(bbox.Left, bbox.Top);
            _topRight = new Point(bbox.Right, bbox.Top);
            _bottomLeft = new Point(bbox.Left, bbox.Bottom);
            _bottomRight = new Point(bbox.Right, bbox.Bottom);
        }

        public double Height
        {
            get { return Compute.EuclideanDistance(_topLeft, _bottomLeft); }
        }

        public double Width
        {
            get { return Compute.EuclideanDistance(_topLeft, _topRight); }
        }

        public PointF Center
        {
            get
            {
                float x = (_topLeft.X + _bottomRight.X) / 2f;
                float y = (_topLeft.Y + _bottomRight.Y) / 2f;

                return new PointF(x, y);
            }
        }

    }

    [Serializable]
    [DebuggerDisplay("angPosition={_angPosition} relDistance={_relDistance}")]
    class PolarPoint : ICloneable
    {
        private double _angPosition;
        private double _relDistance;
        private Guid _id;

        public PolarPoint()
        {
            _id = Guid.NewGuid();
            _angPosition = 0.0;
            _relDistance = 0.0;
        }

        public PolarPoint(double angPosition, double relDistance)
        {
            _id = Guid.NewGuid();
            _angPosition = angPosition;
            _relDistance = relDistance;
        }

        public PolarPoint(PolarPoint pt)
        {
            PolarPoint temp = (PolarPoint)pt.Clone();
            this._id = temp._id;
            this._angPosition = temp._angPosition;
            this._relDistance = temp._relDistance;
        }

        public object Clone()
        {
            PolarPoint pt = (PolarPoint)this.MemberwiseClone();
            return pt;
        }

        public Guid Id
        {
            get { return _id; }
        }

        public double AngularPosition
        {
            get { return _angPosition; }
            set { _angPosition = value; }
        }

        public double RelativeDistance
        {
            get { return _relDistance; }
            set { _relDistance = value; }
        }
    }
}
