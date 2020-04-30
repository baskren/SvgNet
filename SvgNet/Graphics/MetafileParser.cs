﻿/*
    Copyright © 2003 RiskCare Ltd. All rights reserved.
    Copyright © 2010 SvgNet & SvgGdi Bridge Project. All rights reserved.
    Copyright © 2015-2019 Rafael Teixeira, Mojmír Němeček, Benjamin Peterson and Other Contributors

    Original source code licensed with BSD-2-Clause spirit, treat it thus, see accompanied LICENSE for more
*/

using System;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Drawing.Imaging;
using System.IO;

namespace SvgNet.SvgGdi {
    namespace MetafileTools {       // Classes in this namespace were inspired by the code in http://wmf.codeplex.com/
        namespace EmfTools {
            public interface IBinaryRecord {
                void Read(BinaryReader reader);
            }

            public static class BinaryReaderExtensions {
                /// <summary>
                /// Skips excess bytes. Work-around for some WMF files that contain undocumented fields.
                /// </summary>
                /// <param name="reader"></param>
                /// <param name="excess"></param>
                public static void Skip(this BinaryReader reader, int excess) {
                    if (excess > 0) {
                        //Skip unknown bytes
                        reader.BaseStream.Seek(excess, SeekOrigin.Current);
                        //var dummy = reader.ReadBytes(excess);
                    }
                }
            }

            /// <summary>
            /// Implements a EMF META record
            /// </summary>
            public abstract class EmfBinaryRecord : IBinaryRecord {
                /// <summary>
                /// Gets or sets record length
                /// </summary>
                public uint RecordSize {
                    get;
                    set;
                }

                /// <summary>
                /// Gets or sets record type (aka RecordFunction)
                /// </summary>
                public EmfPlusRecordType RecordType {
                    get;
                    set;
                }

                /// <summary>
                /// Reads a record from binary stream. If this method is not overridden it will skip this record and go to next record.
                /// NOTE: When overriding this method remove the base.Read(reader) line from code.
                /// </summary>
                /// <param name="reader"></param>
                public virtual void Read(BinaryReader reader) {
                }

                protected EmfBinaryRecord() {
                }
            }

            [Serializable]
            public sealed class EmfException : Exception {
                public EmfException(string message)
                    : base(message) {
                }

                public EmfException() : base() {
                }

                public EmfException(string message, Exception innerException) : base(message, innerException) {
                }
            }

            /// <summary>
            /// Low-level EMF parser
            /// </summary>
            public class EmfReader : IDisposable {
                public EmfReader(Stream stream) {
                    this._stream = stream;
                    _reader = new BinaryReader(stream);
                }

                public bool IsEndOfFile => _stream.Length == _stream.Position;

                public void Dispose() {
                    if (_reader != null) {
                        _reader.Close();
                        _reader = null;
                    }
                }

                public IBinaryRecord Read() {
                    long begin = _reader.BaseStream.Position;

                    var rt = (EmfPlusRecordType)_reader.ReadUInt32();
                    var recordSize = _reader.ReadUInt32();

                    EmfBinaryRecord record = new EmfUnknownRecord {
                        RecordType = rt,
                        RecordSize = recordSize
                    };
                    record.Read(_reader);

                    long end = _reader.BaseStream.Position;
                    long rlen = end - begin; //Read length
                    long excess = recordSize - rlen;
                    if (excess > 0) {
                        //Oops, reader did not read whole record?!
                        _reader.Skip((int)excess);
                    }

                    return record;
                }

                private readonly Stream _stream;
                private BinaryReader _reader;
            }

            public class EmfUnknownRecord : EmfBinaryRecord {
                public byte[] Data {
                    get;
                    set;
                }

                public override void Read(BinaryReader reader) {
                    var length = (int)base.RecordSize - sizeof(uint) - sizeof(uint);
                    if (length > 0) {
                        Data = reader.ReadBytes(length);
                    } else {
                        Data = _emptyData;
                    }
                }

                private static readonly byte[] _emptyData = new byte[0];
            }
        }

        public delegate void DrawLineDelegate(PointF[] points);

        public delegate void FillPolygonDelegate(PointF[] points, Brush brush);

        public class MetafileParser {
            public enum EmfBrushStyle {
                BS_SOLID = 0x0000,
                BS_NULL = 0x0001,
                BS_HATCHED = 0x0002,
                BS_PATTERN = 0x0003,
                BS_INDEXED = 0x0004,
                BS_DIBPATTERN = 0x0005,
                BS_DIBPATTERNPT = 0x0006,
                BS_PATTERN8X8 = 0x0007,
                BS_DIBPATTERN8X8 = 0x0008,
                BS_MONOPATTERN = 0x0009
            }

            /// <summary>
            /// https://msdn.microsoft.com/en-us/library/cc231191 without the 0x80000000 bit
            /// </summary>
            public enum EmfStockObject {
                WHITE_BRUSH = 0x00000000,
                LTGRAY_BRUSH = 0x00000001,
                GRAY_BRUSH = 0x00000002,
                DKGRAY_BRUSH = 0x00000003,
                BLACK_BRUSH = 0x00000004,
                NULL_BRUSH = 0x00000005,
                WHITE_PEN = 0x00000006,
                BLACK_PEN = 0x00000007,
                NULL_PEN = 0x00000008,
                OEM_FIXED_FONT = 0x0000000A,
                ANSI_FIXED_FONT = 0x0000000B,
                ANSI_VAR_FONT = 0x0000000C,
                SYSTEM_FONT = 0x0000000D,
                DEVICE_DEFAULT_FONT = 0x0000000E,
                DEFAULT_PALETTE = 0x0000000F,
                SYSTEM_FIXED_FONT = 0x00000010,
                DEFAULT_GUI_FONT = 0x00000011,
                DC_BRUSH = 0x00000012,
                DC_PEN = 0x00000013,

                MinValue = WHITE_BRUSH,
                MaxValue = DC_PEN
            }

            public enum EmfTransformMode {
                MWT_IDENTITY = 1,
                MWT_LEFTMULTIPLY = 2,
                MWT_RIGHTMULTIPLY = 3
            }

            public void EnumerateMetafile(Stream emf, float unitSize, PointF destination, DrawLineDelegate drawLine, FillPolygonDelegate fillPolygon) {
                _transform = new SVG.Gaps.Matrix();
                _drawLine = drawLine;
                _fillPolygon = fillPolygon;
                _zero = destination;
                _lineBuffer = new LineBuffer(unitSize);
                _objects = new Dictionary<uint, ObjectHandle>();
                _brush = null;

                using (var reader = new EmfTools.EmfReader(emf)) {
                    while (!reader.IsEndOfFile) if (reader.Read() is EmfTools.EmfUnknownRecord record) Enumerate(record);
                }

                CommitLine();
            }

            private void Enumerate(EmfTools.EmfUnknownRecord record) {
                switch (record.RecordType) {
                    case EmfPlusRecordType.EmfHeader:
                    case EmfPlusRecordType.EmfEof:
                    case EmfPlusRecordType.EmfSaveDC:
                    case EmfPlusRecordType.EmfDeleteObject:
                    case EmfPlusRecordType.EmfExtCreatePen:
                    case EmfPlusRecordType.EmfCreatePen:
                    case EmfPlusRecordType.EmfRestoreDC:
                    case EmfPlusRecordType.EmfSetIcmMode:
                    case EmfPlusRecordType.EmfSetMiterLimit:
                    case EmfPlusRecordType.EmfSetPolyFillMode:
                        // Harmless records with no significant side-effects on the shape of the drawn outline
                        break;

                    case EmfPlusRecordType.EmfSelectObject:
                        ProcessSelectObject(record.Data);
                        break;

                    case EmfPlusRecordType.EmfCreateBrushIndirect:
                        ProcessCreateBrushIndirect(record.Data);
                        break;

                    case EmfPlusRecordType.EmfBeginPath:
                        ProcessBeginPath(record.Data);
                        break;

                    case EmfPlusRecordType.EmfEndPath:
                        // TODO:
                        break;

                    case EmfPlusRecordType.EmfStrokeAndFillPath:
                        ProcessStrokeAndFillPath(record.Data);
                        break;

                    case EmfPlusRecordType.EmfMoveToEx:
                        ProcessMoveToEx(record.Data);
                        break;

                    case EmfPlusRecordType.EmfModifyWorldTransform:
                        ProcessModifyWorldTransform(record.Data);
                        break;

                    case EmfPlusRecordType.EmfPolygon16:
                        ProcessPolygon16(record.Data);
                        break;

                    case EmfPlusRecordType.EmfPolyPolygon16:
                        ProcessPolyPolygon16(record.Data);
                        break;

                    case EmfPlusRecordType.EmfPolyline16:
                        ProcessPolyline16(record.Data);
                        break;

                    case EmfPlusRecordType.EmfPolylineTo16:
                        ProcessPolylineTo16(record.Data);
                        break;

                    case EmfPlusRecordType.EmfCloseFigure:
                        ProcessCloseFigure(record.Data);
                        break;

                    case EmfPlusRecordType.EmfPolyBezierTo16:
                        ProcessPolyBezierTo16(record.Data);
                        break;

                    default:
                        throw new NotImplementedException();
                }
            }

            private const uint _stockObjectMaxCode = 0x80000000 + (uint)EmfStockObject.MaxValue;
            private const uint _stockObjectMinCode = 0x80000000 + (uint)EmfStockObject.MinValue;
            private Brush _brush;
            private PointF _curveOrigin;
            private DrawLineDelegate _drawLine;
            private FillPolygonDelegate _fillPolygon;
            private LineBuffer _lineBuffer;
            private PointF _moveTo;
            private Dictionary<uint, ObjectHandle> _objects;
            private SVG.Gaps.Matrix _transform;
            private PointF _zero;

            private void CommitLine() {
                if (_lineBuffer.IsEmpty)
                    return;

                var linePoints = _lineBuffer.GetPoints();

                _lineBuffer.Clear();

                if (linePoints == null)
                    return;

                for (var i = 0; i < linePoints.Length; i++) {
                    linePoints[i].X += _zero.X;
                    linePoints[i].Y += _zero.Y;
                }
                _drawLine(linePoints);
            }

            private void DrawLine(PointF[] points, int offset, int count) {
                if (!_lineBuffer.CanAdd(points, offset)) {
                    CommitLine();
                }

                _lineBuffer.Add(points, offset, count);
            }

            private void FillPolygon(PointF[] linePoints, Brush fillBrush) {
                if (linePoints == null || fillBrush == null)
                    return;

                for (var i = 0; i < linePoints.Length; i++) {
                    linePoints[i].X += _zero.X;
                    linePoints[i].Y += _zero.Y;
                }
                _fillPolygon(linePoints, fillBrush);
            }

            private void InternalProcessPolyline16(uint numberOfPolygons, uint totalNumberOfPoints, int[] numberOfPoints, BinaryReader reader) {
                var points = new PointF[totalNumberOfPoints];

                for (var j = 0; j < points.Length; j++) {
                    points[j].X = reader.ReadInt16();
                    points[j].Y = reader.ReadInt16();
                }

                _transform.TransformPoints(points);

                var offset = 0;
                for (var i = 0; i < numberOfPolygons; i++) {
                    DrawLine(points, offset, numberOfPoints[i]);
                    offset += numberOfPoints[i];
                }
            }

            private void InternalSelectObject(EmfStockObject stockObject) {
                switch (stockObject) {
                    case EmfStockObject.BLACK_BRUSH:
                        _brush = new SolidBrush(Color.Black);
                        break;

                    case EmfStockObject.DC_BRUSH:
                        throw new NotImplementedException();

                    case EmfStockObject.DKGRAY_BRUSH:
                        _brush = new SolidBrush(Color.DarkGray);
                        break;

                    case EmfStockObject.GRAY_BRUSH:
                        _brush = new SolidBrush(Color.Gray);
                        break;

                    case EmfStockObject.LTGRAY_BRUSH:
                        _brush = new SolidBrush(Color.LightGray);
                        break;

                    case EmfStockObject.NULL_BRUSH:
                        _brush = null;
                        break;

                    case EmfStockObject.WHITE_BRUSH:
                        _brush = new SolidBrush(Color.White);
                        break;
                }
            }

            private void ProcessBeginPath(byte[] recordData) {
                MemoryStream _ms = null;
                BinaryReader _br = null;
                try {
                    _ms = new MemoryStream(recordData);
                    _br = new BinaryReader(_ms);

                    System.Diagnostics.Debug.Assert(_ms.Position == _ms.Length);

                    // Clear the line buffer so that it can record the path
                    CommitLine();
                } finally {
                    _br?.Close();
                    _ms?.Dispose();
                }
            }

            private void ProcessCloseFigure(byte[] recordData) {
                MemoryStream _ms = null;
                BinaryReader _br = null;
                try {
                    _ms = new MemoryStream(recordData);
                    _br = new BinaryReader(_ms);

                    var points = new PointF[2];

                    points[0].X = _moveTo.X;
                    points[0].Y = _moveTo.Y;
                    points[1].X = _curveOrigin.X;
                    points[1].Y = _curveOrigin.Y;

                    _transform.TransformPoints(points);

                    DrawLine(points, 0, points.Length);

                    System.Diagnostics.Debug.Assert(_ms.Position == _ms.Length);
                } finally {
                    _br?.Close();
                    _ms?.Dispose();
                }
            }

            private void ProcessCreateBrushIndirect(byte[] recordData) {
                MemoryStream _ms = null;
                BinaryReader _br = null;
                try {
                    _ms = new MemoryStream(recordData);
                    _br = new BinaryReader(_ms);

                    var ihBrush = _br.ReadUInt32();

                    // https://msdn.microsoft.com/en-us/library/cc230581.aspx
                    var brushStyle = (EmfBrushStyle)_br.ReadUInt32();
                    var r = _br.ReadByte();
                    var g = _br.ReadByte();
                    var b = _br.ReadByte();
                    var reserved = _br.ReadByte();
                    var brushColor = Color.FromArgb(r, g, b);
                    var brushHatch = _br.ReadUInt32();

                    _objects.Remove(ihBrush);

                    switch (brushStyle) {
                        case EmfBrushStyle.BS_SOLID:
                            _objects.Add(ihBrush, new ObjectHandle(new SolidBrush(brushColor)));
                            break;

                        case EmfBrushStyle.BS_NULL:
                            _objects.Add(ihBrush, new ObjectHandle(EmfStockObject.NULL_BRUSH));
                            break;

                        case EmfBrushStyle.BS_HATCHED:
                            throw new NotImplementedException();

                        default:
                            throw new NotImplementedException();
                    }

                    System.Diagnostics.Debug.Assert(_ms.Position == _ms.Length);
                } finally {
                    _br?.Close();
                    _ms?.Dispose();
                }
            }

            private void ProcessModifyWorldTransform(byte[] recordData) {
                MemoryStream _ms = null;
                BinaryReader _br = null;
                try {
                    _ms = new MemoryStream(recordData);
                    _br = new BinaryReader(_ms);

                    var eM11 = _br.ReadSingle();
                    var eM12 = _br.ReadSingle();
                    var eM21 = _br.ReadSingle();
                    var eM22 = _br.ReadSingle();
                    var eDx = _br.ReadSingle();
                    var eDy = _br.ReadSingle();
                    var iMode = (EmfTransformMode)_br.ReadInt32();

                    using (var matrix = new SVG.Gaps.Matrix(eM11, eM12, eM21, eM22, eDx, eDy)) {
                        switch (iMode) {
                            case EmfTransformMode.MWT_IDENTITY:
                                _transform = new SVG.Gaps.Matrix();
                                break;

                            case EmfTransformMode.MWT_LEFTMULTIPLY:
                                _transform.Multiply(matrix, MatrixOrder.Append /* TODO: is it the correct order? */);
                                break;

                            case EmfTransformMode.MWT_RIGHTMULTIPLY:
                                _transform.Multiply(matrix, MatrixOrder.Prepend /* TODO: is it the correct order? */);
                                break;

                            default:
                                throw new NotImplementedException();
                        }
                    }

                    System.Diagnostics.Debug.Assert(_ms.Position == _ms.Length);
                } finally {
                    _br?.Close();
                    _ms?.Dispose();
                }
            }

            private void ProcessMoveToEx(byte[] recordData) {
                MemoryStream _ms = null;
                BinaryReader _br = null;
                try {
                    _ms = new MemoryStream(recordData);
                    _br = new BinaryReader(_ms);

                    _moveTo = new PointF {
                        X = _br.ReadInt32(),
                        Y = _br.ReadInt32()
                    };

                    _curveOrigin = _moveTo;

                    System.Diagnostics.Debug.Assert(_ms.Position == _ms.Length);
                } finally {
                    _br?.Close();
                    _ms?.Dispose();
                }
            }

            private void ProcessPolyBezierTo16(byte[] recordData) {
                MemoryStream _ms = null;
                BinaryReader _br = null;
                try {
                    _ms = new MemoryStream(recordData);
                    _br = new BinaryReader(_ms);

                    _br.ReadBytes(16 /* Bounds */);

                    var totalNumberOfPoints = _br.ReadUInt32();

                    var originalPoints = new PointF[totalNumberOfPoints];

                    for (var j = 0; j < originalPoints.Length; j++) {
                        originalPoints[j].X = _br.ReadInt16();
                        originalPoints[j].Y = _br.ReadInt16();
                    }

                    const int PointsPerCurve = 3;

                    var numberOfCurves = totalNumberOfPoints / PointsPerCurve;

                    var points = new PointF[1 + numberOfCurves];

                    // Clone _moveTo cursor
                    points[0].X = _moveTo.X;
                    points[0].Y = _moveTo.Y;

                    for (var j = 1; j < points.Length; j++) {
                        // Every curve is defined by 3 points. The first two are the Bezier curve's control points.
                        // The 3rd is the endpoint. This is the point we'll use (only)
                        points[j] = originalPoints[((j - 1) * PointsPerCurve) + (PointsPerCurve - 1)];
                    }

                    // Clone last point to the current _moveTo cursor
                    _moveTo = new PointF(points[points.Length - 1].X, points[points.Length - 1].Y);

                    _transform.TransformPoints(points);

                    DrawLine(points, 0, points.Length);

                    System.Diagnostics.Debug.Assert(_ms.Position == _ms.Length);
                } finally {
                    _br?.Close();
                    _ms?.Dispose();
                }
            }

            private void ProcessPolygon16(byte[] recordData) {
                MemoryStream _ms = null;
                BinaryReader _br = null;
                try {
                    _ms = new MemoryStream(recordData);
                    _br = new BinaryReader(_ms);

                    _br.ReadBytes(16 /* Bounds */);

                    var totalNumberOfPoints = _br.ReadUInt32();

                    var numberOfPoints = new int[1];
                    numberOfPoints[0] = (int)totalNumberOfPoints;

                    InternalProcessPolyline16(1, totalNumberOfPoints, numberOfPoints, _br);

                    System.Diagnostics.Debug.Assert(_ms.Position == _ms.Length);
                } finally {
                    _br?.Close();
                    _ms?.Dispose();
                }
            }

            private void ProcessPolyline16(byte[] recordData) {
                MemoryStream _ms = null;
                BinaryReader _br = null;
                try {
                    _ms = new MemoryStream(recordData);
                    _br = new BinaryReader(_ms);

                    _br.ReadBytes(16 /* Bounds */);

                    const uint numberOfPolygons = 1u;
                    var totalNumberOfPoints = _br.ReadUInt32();

                    var numberOfPoints = new int[numberOfPolygons];
                    numberOfPoints[0] = (int)totalNumberOfPoints;

                    InternalProcessPolyline16(numberOfPolygons, totalNumberOfPoints, numberOfPoints, _br);

                    System.Diagnostics.Debug.Assert(_ms.Position == _ms.Length);
                } finally {
                    _br?.Close();
                    _ms?.Dispose();
                }
            }

            private void ProcessPolylineTo16(byte[] recordData) {
                MemoryStream _ms = null;
                BinaryReader _br = null;
                try {
                    _ms = new MemoryStream(recordData);
                    _br = new BinaryReader(_ms);

                    _br.ReadBytes(16 /* Bounds */);

                    var totalNumberOfPoints = _br.ReadUInt32();

                    var points = new PointF[1 + totalNumberOfPoints];

                    // Clone _moveTo cursor
                    points[0].X = _moveTo.X;
                    points[0].Y = _moveTo.Y;

                    for (var j = 1; j < points.Length; j++) {
                        points[j].X = _br.ReadInt16();
                        points[j].Y = _br.ReadInt16();
                    }

                    // Clone last point to the current _moveTo cursor
                    _moveTo = new PointF(points[points.Length - 1].X, points[points.Length - 1].Y);

                    _transform.TransformPoints(points);

                    DrawLine(points, 0, points.Length);

                    System.Diagnostics.Debug.Assert(_ms.Position == _ms.Length);
                } finally {
                    _br?.Close();
                    _ms?.Dispose();
                }
            }

            private void ProcessPolyPolygon16(byte[] recordData) {
                MemoryStream _ms = null;
                BinaryReader _br = null;
                try {
                    _ms = new MemoryStream(recordData);
                    _br = new BinaryReader(_ms);

                    _br.ReadBytes(16 /* Bounds */);

                    var numberOfPolygons = _br.ReadUInt32();
                    var totalNumberOfPoints = _br.ReadUInt32();

                    var numberOfPoints = new int[numberOfPolygons];

                    for (var i = 0; i < numberOfPolygons; i++)
                        numberOfPoints[i] = (int)_br.ReadUInt32();

                    InternalProcessPolyline16(numberOfPolygons, totalNumberOfPoints, numberOfPoints, _br);

                    System.Diagnostics.Debug.Assert(_ms.Position == _ms.Length);
                } finally {
                    _br?.Close();
                    _ms?.Dispose();
                }
            }

            private void ProcessSelectObject(byte[] recordData) {
                MemoryStream _ms = null;
                BinaryReader _br = null;
                try {
                    _ms = new MemoryStream(recordData);
                    _br = new BinaryReader(_ms);

                    var ihObject = _br.ReadUInt32();

                    if (ihObject >= _stockObjectMinCode && ihObject <= _stockObjectMaxCode) {
                        var stockObject = (EmfStockObject)(ihObject - _stockObjectMinCode + (int)EmfStockObject.MinValue);
                        InternalSelectObject(stockObject);
                    } else {
                        if (_objects.TryGetValue(ihObject, out ObjectHandle objectHandle)) {
                            if (objectHandle.IsStockObject) {
                                InternalSelectObject(objectHandle.GetStockObject());
                            } else if (objectHandle.IsBrush) {
                                _brush = objectHandle.GetBrush();
                            }
                        }
                    }

                    System.Diagnostics.Debug.Assert(_ms.Position == _ms.Length);
                } finally {
                    _br?.Close();
                    _ms?.Dispose();
                }
            }

            private void ProcessStrokeAndFillPath(byte[] recordData) {
                MemoryStream _ms = null;
                BinaryReader _br = null;
                try {
                    _ms = new MemoryStream(recordData);
                    _br = new BinaryReader(_ms);

                    _br.ReadBytes(16 /* Bounds */);

                    System.Diagnostics.Debug.Assert(_ms.Position == _ms.Length);

                    FillPolygon(_lineBuffer.GetPoints(), _brush);
                } finally {
                    _br?.Close();
                    _ms?.Dispose();
                }
            }

            private class LineBuffer {
                public LineBuffer(float unitSize) {
                    _points = new List<NormalizedPoint>();
                    _normalizedPoints = new List<NormalizedPoint>();
                    _visualPoints = new List<VisualPoint>();
                    _epsilonSquare = _unitSizeEpsilon * unitSize * _unitSizeEpsilon * unitSize;
                }

                public bool IsEmpty => _points.Count == 0;

                public void Add(PointF[] points, int offset, int count) {
                    if (IsEmpty) {
                        MakeRoom(count);
                        for (var i = 0; i < count; i++)
                            _points.Add(Add(points[offset + i]));
                    } else {
                        if (!IsVisuallyIdentical(GetLastPoint(), points[offset]))
                            throw new ArgumentOutOfRangeException();

                        MakeRoom(count - 1);

                        Add(points[offset]);

                        for (var i = 1; i < count; i++)
                            _points.Add(Add(points[offset + i]));
                    }
                }

                public bool CanAdd(PointF[] points, int offset) {
                    if (IsEmpty)
                        return true;

                    if (IsVisuallyIdentical(GetLastPoint(), points[offset]))
                        return true;

                    return false;
                }

                public void Clear() => _points.Clear();

                public PointF[] GetPoints() {
                    var points = new List<NormalizedPoint> {
                        _points[0]
                    };
                    for (var i = 1; i < _points.Count; i++) {
                        if (!IsVisuallyIdentical(points[points.Count - 1], _points[i])) {
                            points.Add(_points[i]);
                        }
                    }

                    if (points.Count <= 1)
                        return null;

                    var result = new List<PointF>();

                    for (var i = 0; i < points.Count; i++) {
                        var visualPoint = _visualPoints[points[i].VisualIndex];
                        if (!visualPoint.IsLocked) {
                            // Calculate the visual point's appearance as "the middle" of all points
                            double sumX = 0;
                            double sumY = 0;
                            for (int j = 0, siblingCount = visualPoint.Weight; siblingCount > 0; j++) {
                                if (_normalizedPoints[j].VisualIndex == visualPoint.VisualIndex) {
                                    sumX += _normalizedPoints[j].Point.X;
                                    sumY += _normalizedPoints[j].Point.Y;
                                    siblingCount--;
                                }
                            }

                            visualPoint.Point = new PointF((float)(sumX / visualPoint.Weight), (float)(sumY / visualPoint.Weight));
                            visualPoint.IsLocked = true;
                        }
                        result.Add(visualPoint.Point);
                    }

                    return result.ToArray();
                }

                private const float _unitSizeEpsilon = 2.0f;
                private readonly float _epsilonSquare;
                private readonly List<NormalizedPoint> _normalizedPoints;
                private readonly List<NormalizedPoint> _points;
                private readonly List<VisualPoint> _visualPoints;

                private static bool IsVisuallyIdentical(NormalizedPoint a, NormalizedPoint b) => a.VisualIndex == b.VisualIndex;

                private NormalizedPoint Add(PointF point) {
                    NormalizedPoint result;
                    VisualPoint visualPoint;

                    for (var i = _normalizedPoints.Count - 1; i >= 0; i--) {
                        if (IsVisuallyIdentical(_normalizedPoints[i].Point, point)) {
                            visualPoint = _visualPoints[_normalizedPoints[i].VisualIndex];
                            visualPoint.Weight++;
                            result = new NormalizedPoint { Point = point, VisualIndex = visualPoint.VisualIndex };
                            _normalizedPoints.Add(result);
                            return result;
                        }
                    }

                    visualPoint = new VisualPoint { IsLocked = false, VisualIndex = _visualPoints.Count, Weight = 1 };
                    _visualPoints.Add(visualPoint);

                    result = new NormalizedPoint { Point = point, VisualIndex = visualPoint.VisualIndex };
                    _normalizedPoints.Add(result);
                    return result;
                }

                private NormalizedPoint GetLastPoint() => _points[_points.Count - 1];

                // TODO: TUNE: what's the correct value? Shoud it be based on the matafile's DPI?
                private bool IsVisuallyIdentical(NormalizedPoint a, PointF b) {
                    for (var i = _normalizedPoints.Count - 1; i >= 0; i--)
                        if (_normalizedPoints[i].VisualIndex == a.VisualIndex && IsVisuallyIdentical(_normalizedPoints[i].Point, b)) return true;
                    return false;
                }

                private bool IsVisuallyIdentical(PointF a, PointF b) {
                    var dx = a.X - b.X;
                    var dy = a.Y - b.Y;
                    return ((dx * dx) + (dy * dy)) <= _epsilonSquare;
                }

                private void MakeRoom(int count) {
                    if (_points.Capacity < _points.Count + count)
                        _points.Capacity = _points.Count + count;
                }
            }

            private class NormalizedPoint {
                public PointF Point;
                public int VisualIndex;
            }

            private class ObjectHandle {
                public ObjectHandle(EmfStockObject stockObject) => _stockObject = stockObject;

                public ObjectHandle(Brush brush) => _brush = brush;

                public bool IsBrush => _brush != null;

                public bool IsStockObject => _stockObject != null;

                public Brush GetBrush() => _brush;

                public EmfStockObject GetStockObject() => _stockObject.Value;

                private readonly Brush _brush;
                private readonly EmfStockObject? _stockObject;
            }

            private class VisualPoint {
                public bool IsLocked;
                public PointF Point;
                public int VisualIndex;
                public int Weight;
            }
        }
    }
}
