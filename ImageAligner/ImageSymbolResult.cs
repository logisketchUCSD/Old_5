using System;
using System.Collections.Generic;
using System.Text;

namespace ImageAligner
{
    [Serializable]
    public class ImageSymbolResult
    {
        /// <summary>
        /// Scaling factor to keep scores in a nice range
        /// </summary>
        public static double HAUS_SCALE = 1.0 / 48.0;

        /// <summary>
        /// Scaling factor to keep scores in a nice range
        /// </summary>
        public static double MODHAUS_SCALE = 2.0 / 48.0;

        #region Member Variables

        ImageSymbol _UnknownSymbol;

        ImageSymbol _TemplateSymbol;

        double _HausdorffDistance;

        double _ModifiedHausdorffDistance;

        double _TanimotoCoefficient;

        double _YuleCoefficient;

        double _Score = -1.0;

        #endregion

        public ImageSymbolResult()
        {
        }

        public ImageSymbolResult(ImageSymbol unknown, ImageSymbol template, double haus, double modHaus, double tanimoto, double yule)
        {
            _UnknownSymbol = unknown;
            _TemplateSymbol = template;
            _HausdorffDistance = haus;
            _ModifiedHausdorffDistance = modHaus;
            _TanimotoCoefficient = tanimoto;
            _YuleCoefficient = yule;

            CalculateScore();
        }

        private void CalculateScore()
        {
            _Score = 1.0 / (1.0 + _HausdorffDistance * HAUS_SCALE + _ModifiedHausdorffDistance * MODHAUS_SCALE);
        }

        public double Score
        {
            get { return _Score; }
        }

        public string Name
        {
            get { return _TemplateSymbol.Name; }
        }

        public double Hausdorff
        {
            get { return _HausdorffDistance; }
        }

        public double ModifiedHausdorff
        {
            get { return _ModifiedHausdorffDistance; }
        }
    }
}
