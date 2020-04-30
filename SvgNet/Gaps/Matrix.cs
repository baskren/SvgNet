using System;
using System.Drawing;
using System.Drawing.Drawing2D;

namespace SVG.Gaps
{
	public sealed class Matrix : MarshalByRefObject, IDisposable {

		static float[] _identity = new float[] { 1, 0, 0, 1, 0, 0 };
		MatrixTypes _type = MatrixTypes.Identity;
		int _padding = 0;

		float[] _elements = new float[] { 1, 0, 0, 1, 0, 0 };
		public float[] Elements => _elements;

		float _m11 {
			get => _elements[0];
			set => _elements[0] = value;
        }

		float _m12 {
			get => _elements[1];
			set => _elements[1] = value;
		}

		float _m21 {
			get => _elements[2];
			set => _elements[2] = value;
		}

		float _m22 {
			get => _elements[3];
			set => _elements[3] = value;
		}

		float _dx {
			get => _elements[4];
			set => _elements[4] = value;
		}

		float _dy {
			get => _elements[5];
			set => _elements[5] = value;
		}


		public bool IsIdentity =>
			Elements[0] == 1 && Elements[1] == 0 &&
			Elements[2] == 0 && Elements[3] == 1 &&
			Elements[4] == 0 && Elements[5] == 0;

		public bool IsInvertible {
			get {
				if (_m12 == 0 && _m21 == 0)
					return _m11 != 0 && _m22 != 0;
				return Math.Abs(Determinant()) >= 1e-5;
			}
		}

		public float OffsetX => Elements[4];

		public float OffsetY => Elements[5];

		public Matrix() {
		}

		public Matrix(Rectangle rect, Point[] pt) {
			if (pt is null)
				throw new ArgumentNullException(nameof(pt));
			if (pt.Length != 3)
				throw new ArgumentException(nameof(pt));

			_m11 = (pt[1].X - pt[0].X) / rect.Width;
			_m21 = (pt[2].X - pt[0].X) / rect.Height;
			_dx = pt[0].X - _m11 * rect.X - _m21 * rect.Y;
			_m12 = (pt[1].Y - pt[0].Y) / rect.Width;
			_m22 = (pt[2].Y - pt[0].Y) / rect.Height;
			_dy = pt[0].Y -_m12 * rect.X - _m22 * rect.Y;
		}

		public Matrix(RectangleF rect, PointF[] pt) {
			if (pt is null)
				throw new ArgumentNullException(nameof(pt));
			if (pt.Length != 3)
				throw new ArgumentException(nameof(pt));

			_m11 = (pt[1].X - pt[0].X) / rect.Width;
			_m21 = (pt[2].X - pt[0].X) / rect.Height;
			_dx = pt[0].X - _m11 * rect.X - _m21 * rect.Y;
			_m12 = (pt[1].Y - pt[0].Y) / rect.Width;
			_m22 = (pt[2].Y - pt[0].Y) / rect.Height;
			_dy = pt[0].Y - _m12 * rect.X - _m22 * rect.Y;
		}

		public Matrix(float m11, float m12, float m21, float m22, float dx, float dy)
        {
			_m11 = m11;
			_m12 = m12;
			_m21 = m21;
			_m22 = m22;
			_dx = dx;
			_dy = dy;
			_type = MatrixTypes.Unknown;
			DeriveType();
		}

		public Matrix Clone()
			=> new Matrix(Elements[0], Elements[1], Elements[2], Elements[3], Elements[4], Elements[5]);

		public void Dispose() {
		}

		public override bool Equals(object obj) {
			if (obj is Matrix other) {
				return
					Elements[0] == other.Elements[0] &&
					Elements[1] == other.Elements[1] &&
					Elements[2] == other.Elements[2] &&
					Elements[3] == other.Elements[3] &&
					Elements[4] == other.Elements[4] &&
					Elements[5] == other.Elements[5];
			}
			return false;
		}

		public void Invert() {
			if (!IsInvertible)
				return;
            if (_m12 == 0 && _m21 == 0) {
				_dx = -_dx / _m11;
				_dy = -_dy / _m22;
				_m11 = 1 / _m11;
				_m22 = 1 / _m22;
            }
            else {
				var det = Determinant();
				var copy = Clone();
				_m11 = copy._m22 / det;
				_m12 = -copy._m12 / det;
				_m21 = -copy._m21 / det;
				_m22 = copy._m11 / det;
				_dx = (copy._m21 * copy._dy - copy._m22 * copy._dx) / det;
				_dy = (copy._m11 * copy._dy - copy._m12 * copy._dx) / det;
            }
		}

		public void Multiply(Matrix matrix) 
			=> Multiply(matrix, MatrixOrder.Prepend);
		

		public void Multiply(Matrix matrix, MatrixOrder order) {
			if (matrix is null)
				return;

			if (order == MatrixOrder.Append)
				_elements = Mul(this._elements, matrix._elements);
			else
				_elements = Mul(matrix._elements, this._elements);
		}


        static float[] Mul(float[] a, float[] b) {
			var r = new float[6];

			//r[0] = a._m11 * b._m11 + a._m12 * b._m21;
			r[0] = a[0] * b[0] + a[1] * b[2];
			//r[1] = a._m11 * b._m12 + a._m12 * b._m22;
			r[1] = a[0] * b[1] + a[1] * b[3];
			//r[2] = a._m21 * b._m11 + a._m22 * b._m21;
			r[2] = a[2] * b[0] + a[3] * b[2];
			//r[3] = a._m21 * b._m12 + a._m22 * b._m22;
			r[3] = a[2] * b[1] + a[3] + b[3];
			//r[4] = a._dx * b._m11 + a._dy * b._m21 + b._dx;
			r[4] = a[4] * b[0] + a[5] * b[2] + b[4];
			//r[5] = a._dx * b._m12 + a._dy * b._m22 + b._dy;
			r[5] = a[4] * b[1] + a[5] * b[3] + b[5];
			return r;
        }

		public void Reset() {
			_elements = new float[] { 1, 0, 0, 1, 0, 0 };
		}

		public void Rotate(float angle)
            => Rotate(angle, MatrixOrder.Prepend);

		public void Rotate(float angle, MatrixOrder order) {
			float cos_theta, sin_theta;
			float[] rotate = new float[6];

			angle *= (float)(Math.PI / 180.0);
			cos_theta = (float)Math.Cos(angle);
			sin_theta = (float)Math.Sin(angle);
            rotate[0] = cos_theta;
			rotate[1] = sin_theta;
			rotate[2] = -sin_theta;
			rotate[3] = cos_theta;
			rotate[4] = 0.0f;
			rotate[5] = 0.0f;

            if (order == MatrixOrder.Append)
				_elements = Mul(_elements, rotate);
            else
				_elements = Mul(rotate, _elements);
		}

		public void RotateAt(float angle, PointF point) 
			=> RotateAt(angle, point, MatrixOrder.Prepend);
		

		public void RotateAt(float angle, PointF point, MatrixOrder order) {
			angle *= (float)(Math.PI / 180.0);  // degrees to radians
			var cos = (float)Math.Cos(angle);
			var sin = (float)Math.Sin(angle);
			var e4 = -point.X * cos + point.Y * sin + point.X;
			var e5 = -point.X * sin - point.Y * cos + point.Y;
			var m = new float[] { _m11, _m12, _m21, _m22, _dx, _dy };

			if (order == MatrixOrder.Prepend) {
				_m11 = cos * m[0] + sin * m[2];
				_m12 = cos * m[1] + sin * m[3];
				_m12 = -sin * m[0] + cos * m[2];
				_m22 = -sin * m[1] + cos * m[3];
				_dx = e4 * m[0] + e5 * m[2] + m[4];
				_dy = e4* m[1] +e5 * m[3] + m[5];
			} else {
				_m11 = m[0] * cos + m[1] * -sin;
				_m12 = m[0] * sin + m[1] * cos;
				_m21 = m[2] * cos + m[3] * -sin;
				_m22 = m[2] * sin + m[3] * cos;
				_dx = m[4] * cos + m[5] * -sin + e4;
				_dy = m[4] * sin + m[5] * cos + e5;
			}
		}

		public void Scale(float scaleX, float scaleY) 
			=> Scale(scaleX, scaleY, MatrixOrder.Prepend);
		

		public void Scale(float scaleX, float scaleY, MatrixOrder order) {
			var scale = new float[] { scaleX, 0, 0, scaleY, 0, 0 };

			if (order == MatrixOrder.Append)
				_elements = Mul(_elements, scale);
			else
				_elements = Mul(scale, _elements);
		}

		public void Shear(float shearX, float shearY) 
			=> Shear(shearX, shearY, MatrixOrder.Prepend);
		

		public void Shear(float shearX, float shearY, MatrixOrder order) {
            var shear = new float[] { 1, shearY, shearX, 1, 0, 0  };

			if (order == MatrixOrder.Append)
				_elements = Mul(_elements, shear);
			else
				_elements = Mul(shear, _elements);
		}

		public void TransformPoints(PointF[] pts) {
			if (pts is null || pts.Length < 1)
				throw new ArgumentNullException(nameof(pts));

			float x, y;
            for(int i=0;i<pts.Length;i++) {
				var pt = pts[i];
				x = pt.X;
				y = pt.Y;
				pt.X = x * _m11 + y * _m21 + _dx;
				pt.Y = x * _m12 + y * _m22 + _dy;
            }
		}

		public void TransformPoints(Point[] pts) {
			if (pts is null || pts.Length < 1)
				throw new ArgumentNullException(nameof(pts));

			float x, y;
			for (int i = 0; i < pts.Length; i++) {
				var pt = pts[i];
				x = pt.X;
				y = pt.Y;
				pt.X = (int)Math.Round(x * _m11 + y * _m21 + _dx);
				pt.Y = (int)Math.Round(x * _m12 + y * _m22 + _dy);
			}
		}

		public void TransformVectors(PointF[] pts) {
			if (pts is null || pts.Length < 1)
				throw new ArgumentNullException(nameof(pts));

			float x, y;
			for (int i = 0; i < pts.Length; i++) {
				var pt = pts[i];
				x = pt.X;
				y = pt.Y;
				pt.X = x * _m11 + y * _m21;
				pt.Y = x * _m12 + y * _m22;
			}
		}

		public void TransformVectors(Point[] pts) {
			if (pts is null || pts.Length < 1)
				throw new ArgumentNullException(nameof(pts));

			float x, y;
			for (int i = 0; i < pts.Length; i++) {
				var pt = pts[i];
				x = pt.X;
				y = pt.Y;
				pt.X = (int)Math.Round(x * _m11 + y * _m21);
				pt.Y = (int)Math.Round(x * _m12 + y * _m22);
			}
		}

		public void Translate(float offsetX, float offsetY) {
			Translate(offsetX, offsetY, MatrixOrder.Prepend);
		}

		public void Translate(float offsetX, float offsetY, MatrixOrder order) {
			var trans = new float[] { 1,0,0,1,offsetX, offsetY };

			if (order == MatrixOrder.Append)
				_elements = Mul(_elements, trans);
			else
				_elements = Mul(trans, _elements);
		}

		public void VectorTransformPoints(Point[] pts) {
			TransformVectors(pts);
		}

        void DeriveType()
        {
			_type = 0;
            if (!(_m21 ==0 && _m12 == 0))
            {
				_type = MatrixTypes.Unknown;
				return;
            }
			if (!(_m11 == 1 && _m22 == 1))
				_type = MatrixTypes.Scaling;
			if (!(_dx == 0 && _dy == 0))
				_type |= MatrixTypes.Translation;
			if (0 == (_type & (MatrixTypes.Translation | MatrixTypes.Scaling)))
				_type = MatrixTypes.Identity;
        }

        float Determinant() {
			return (float)(_m11 * (double)_m22 - _m12 * (double)_m21);
        }

        public override int GetHashCode() {
            int hashCode = -1180554807;
            hashCode = hashCode * -1521134295 + _m11.GetHashCode();
            hashCode = hashCode * -1521134295 + _m12.GetHashCode();
            hashCode = hashCode * -1521134295 + _m22.GetHashCode();
            hashCode = hashCode * -1521134295 + _dx.GetHashCode();
            hashCode = hashCode * -1521134295 + _dy.GetHashCode();
            return hashCode;
        }
    }
}