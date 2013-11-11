/*
 * File: StrokeMapping.cs
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
    /// 
    /// </summary>
    [DebuggerDisplay("{m_Map.Count} Pairs, Un-Mapped: {m_UnMappedKeys.Count} in unknown, {m_UnMappedValues.Count} in template")]
    [Serializable]
    public class StrokeMapping
    {
        #region Constants

        /// <summary>
        /// Hash value used in combination with keys (strokes in the 'unknown')
        /// </summary>
        const int hashKey = 1234;

        /// <summary>
        /// Hash value used in combination with values (strokes in the 'known' template)
        /// </summary>
        const int hashValue = 4321;

        #endregion


        #region Member Variables

        /// <summary>
        /// Overall map, indexed by the combined hash code of each pair
        /// </summary>
        Dictionary<int, KeyValuePair<Substroke, Substroke>> m_Map;

        /// <summary>
        /// Hash code of an unknown stroke to its pair
        /// </summary>
        Dictionary<int, int> m_Key2Map;

        /// <summary>
        /// Hash code of a template stroke to its pair
        /// </summary>
        Dictionary<int, int> m_Value2Map;

        /// <summary>
        /// Hash code of an unknown stroke to its paired template stroke
        /// </summary>
        Dictionary<int, int> m_Key2Value;

        /// <summary>
        /// Hash code of a template stroke to its paired unknown stroke
        /// </summary>
        Dictionary<int, int> m_Value2Key;

        /// <summary>
        /// List of all hash codes for unknown strokes
        /// </summary>
        List<int> m_Keys;

        /// <summary>
        /// List of all hash codes for template strokes
        /// </summary>
        List<int> m_Values;

        /// <summary>
        /// List of all strokes in the unknown template that have not been paired
        /// </summary>
        List<Substroke> m_UnMappedKeys;

        /// <summary>
        /// List of all the strokes in the 'known' template that have not been paired
        /// </summary>
        List<Substroke> m_UnMappedValues;

        /// <summary>
        /// Lookup table of a stroke to its type
        /// </summary>
        Dictionary<Substroke, string> m_ShapeParts;

        /// <summary>
        /// Hash of complete mapping, used to quickly compare mappings
        /// </summary>
        int m_HashCode;

        #endregion


        #region Constructors

        /// <summary>
        /// Default Constructor, creates an empty map
        /// </summary>
        public StrokeMapping()
        {
            m_Map = new Dictionary<int, KeyValuePair<Substroke, Substroke>>();
            m_Key2Map = new Dictionary<int, int>();
            m_Value2Map = new Dictionary<int, int>();
            m_Key2Value = new Dictionary<int, int>();
            m_Value2Key = new Dictionary<int, int>();
            m_Keys = new List<int>();
            m_Values = new List<int>();
            m_UnMappedKeys = new List<Substroke>();
            m_UnMappedValues = new List<Substroke>();
            m_ShapeParts = new Dictionary<Substroke, string>();
            m_HashCode = GetHashCode();
        }

        /// <summary>
        /// Clones a stroke mapping
        /// </summary>
        /// <param name="map"></param>
        public StrokeMapping(StrokeMapping map)
        {
            m_Map = new Dictionary<int, KeyValuePair<Substroke, Substroke>>(map.m_Map);
            m_Key2Map = new Dictionary<int, int>(map.m_Key2Map);
            m_Value2Map = new Dictionary<int, int>(map.m_Value2Map);
            m_Key2Value = new Dictionary<int, int>(map.m_Key2Value);
            m_Value2Key = new Dictionary<int, int>(map.m_Value2Key);
            m_Keys = new List<int>(map.m_Keys);
            m_Values = new List<int>(map.m_Values);
            m_UnMappedKeys = new List<Substroke>(map.m_UnMappedKeys);
            m_UnMappedValues = new List<Substroke>(map.m_UnMappedValues);
            m_ShapeParts = map.m_ShapeParts;
            m_HashCode = map.m_HashCode;
        }

        #endregion


        #region Functions

        /// <summary>
        /// Adds a new mapping pair to the StrokeMapping
        /// </summary>
        /// <param name="key">Substroke in the 'unknown'</param>
        /// <param name="value">Substroke in the 'known' template</param>
        public void Add(Substroke key, Substroke value)
        {
            int hash1 = key.Id.GetHashCode();
            int hash2 = value.Id.GetHashCode();
            int hash = hash1 ^ hashKey + hash2 ^ hashValue;
            if (m_Map.ContainsKey(hash))
                return;

            try
            {
                m_Map.Add(hash, new KeyValuePair<Substroke, Substroke>(key, value));
                
                m_Keys.Add(hash1);
                m_Values.Add(hash2);

                m_Key2Map.Add(hash1, hash);
                m_Value2Map.Add(hash2, hash);

                m_Key2Value.Add(hash1, hash2);
                m_Value2Key.Add(hash2, hash1);

                m_UnMappedKeys.Remove(key);
                m_UnMappedValues.Remove(value);
                
                m_HashCode = GetHashCode();
            }
            catch (Exception e)
            {
                throw e;
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="unknown"></param>
        /// <param name="template"></param>
        public void DetermineUnmapped(ImageTemplate unknown, ImageTemplate template)
        {
            foreach (Substroke stroke in unknown.Strokes)
            {
                if (!m_ShapeParts.ContainsKey(stroke))
                    m_ShapeParts.Add(stroke, unknown.GetShapePart(stroke));
                if (!ContainsKey(stroke) && !m_UnMappedKeys.Contains(stroke))
                    m_UnMappedKeys.Add(stroke);
            }

            foreach (Substroke stroke in template.Strokes)
            {
                if (!m_ShapeParts.ContainsKey(stroke))
                    m_ShapeParts.Add(stroke, template.GetShapePart(stroke));
                if (!ContainsValue(stroke) && !m_UnMappedValues.Contains(stroke))
                    m_UnMappedValues.Add(stroke);
            }
        }

        public void RemoveKey(Substroke key)
        {
            try
            {
                int hashKey = key.Id.GetHashCode();
                int hashValue = m_Key2Value[hashKey];
                int hashMap = m_Key2Map[hashKey];

                m_Keys.Remove(hashKey);
                m_Key2Map.Remove(hashKey);
                m_Key2Value.Remove(hashKey);

                m_Values.Remove(hashValue); 
                m_Value2Key.Remove(hashValue);
                m_Value2Map.Remove(hashValue);

                m_Map.Remove(hashMap);

                m_HashCode = GetHashCode();
            }
            catch (Exception e)
            {
                throw e;
            }
        }

        public void RemovePair(Substroke key, Substroke value)
        {
            try
            {
                int hashKey = key.Id.GetHashCode();
                int hashValue = value.Id.GetHashCode();
                int hashMap = m_Key2Map[hashKey];

                m_Keys.Remove(hashKey);
                m_Key2Map.Remove(hashKey);
                m_Key2Value.Remove(hashKey);

                m_Values.Remove(hashValue);
                m_Value2Key.Remove(hashValue);
                m_Value2Map.Remove(hashValue);

                m_Map.Remove(hashMap);

                m_HashCode = GetHashCode();
            }
            catch (Exception e)
            {
                throw e;
            }
        }

        public void RemoveValue(Substroke value)
        {
            try
            {
                int hashValue = value.Id.GetHashCode();
                int hashKey = m_Value2Key[hashValue];
                int hashMap = m_Key2Map[hashKey];

                m_Keys.Remove(hashKey);
                m_Key2Map.Remove(hashKey);
                m_Key2Value.Remove(hashKey);

                m_Values.Remove(hashValue);
                m_Value2Key.Remove(hashValue);
                m_Value2Map.Remove(hashValue);

                m_Map.Remove(hashMap);

                m_HashCode = GetHashCode();
            }
            catch (Exception e)
            {
                throw e;
            }
        }

        public bool ContainsKey(Substroke stroke)
        {
            int hash = stroke.Id.GetHashCode();

            return m_Keys.Contains(hash);
        }

        public bool ContainsValue(Substroke stroke)
        {
            int hash = stroke.Id.GetHashCode();

            return m_Values.Contains(hash);
        }

        public bool ContainsMap(Substroke key, Substroke value)
        {
            int hash = key.Id.GetHashCode() ^ hashKey + value.Id.GetHashCode() ^ hashValue;
            
            return m_Map.ContainsKey(hash);
        }

        #endregion


        #region Getters

        public List<Substroke> StrokeKeys
        {
            get
            {
                List<Substroke> strokes = new List<Substroke>();

                foreach (KeyValuePair<Substroke, Substroke> kv in m_Map.Values)
                    strokes.Add(kv.Key);

                return strokes;
            }
        }

        public List<Substroke> StrokeValues
        {
            get
            {
                List<Substroke> strokes = new List<Substroke>();

                foreach (KeyValuePair<Substroke, Substroke> kv in m_Map.Values)
                    strokes.Add(kv.Value);

                return strokes;
            }
        }

        public List<KeyValuePair<Substroke, Substroke>> Pairs
        {
            get { return new List<KeyValuePair<Substroke, Substroke>>(m_Map.Values); }
        }

        public int Count
        {
            get { return m_Map.Count; }
        }

        public List<KeyValuePair<string, string>> MapWithPartNames
        {
            get
            {
                List<KeyValuePair<string, string>> parts = new List<KeyValuePair<string, string>>();

                foreach (KeyValuePair<Substroke, Substroke> pair in m_Map.Values)
                {
                    string key = m_ShapeParts[pair.Key];
                    string value = m_ShapeParts[pair.Value];
                    parts.Add(new KeyValuePair<string, string>(key, value));
                }

                foreach (Substroke stroke in m_UnMappedValues)
                    parts.Add(new KeyValuePair<string, string>("UnMapped", m_ShapeParts[stroke]));

                foreach (Substroke stroke in m_UnMappedKeys)
                    parts.Add(new KeyValuePair<string, string>(m_ShapeParts[stroke], "UnMapped"));

                return parts;
            }
        }

        public List<Substroke> UnMappedStrokeValues
        {
            get { return m_UnMappedValues; }
        }

        #endregion


        #region Overrides

        public override bool Equals(object obj)
        {
            int hash = m_HashCode;

            StrokeMapping other = (StrokeMapping)obj;
            int otherHash = other.m_HashCode;

            return (hash == otherHash);
        }

        public override int GetHashCode()
        {
            int hash = 0;
            foreach (int h in m_Map.Keys)
                hash += h;

            return hash;
        }

        #endregion
    }

    class StrokeLink
    {
        Substroke m_Key;
        Substroke m_Value;
        int m_Hash;

        public StrokeLink(Substroke key, Substroke value)
        {
            m_Key = key;
            m_Value = value;
            m_Hash = key.Id.GetHashCode() ^ "Key".GetHashCode() + value.Id.GetHashCode() ^ "Value".GetHashCode();
        }

        public int Hash
        {
            get { return m_Hash; }
        }

        public Substroke Key
        {
            get { return m_Key; }
        }

        public Substroke Value
        {
            get { return m_Value; }
        }
    }
}
