using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace NihilKri {
	public class KN {
		public static string Test() {
			//Complex a = new Complex(2.482, -2.736);
			//Complex a = Complex.Cis(1, Math.PI / 3.0 * 5.0);
			//Quaternion A = new Quaternion(8, 8, 8, 8), B = new Quaternion(5, 6, 7, 8);
			//return (Quaternion.k*Quaternion.j*Quaternion.i).ToAbi();
			//return (Math.E ^ new Complex(0, Math.PI)).ToString();
			//return (a^(a*Complex.i)).ToString();
			//return A.ToString();
			//return (A * B / A / B).ToString();

			//Quaternion s = new Quaternion(0, 0, 0, 0), q = A; int nf = 1;
			//for(int n=0 ; n < 20 ; n++) {
			//	s = s + (q / nf);
			//	nf *= n + 1; q = q * A;
			//}
			//return s.ToString();

			//Matrix2 AA = new Matrix2(4, 2, 1, 2, 3, 4, 5, 6, 7, 8);
			//Matrix2 BB = new Matrix2(2, 4, 8, 7, 6, 5, 4, 3, 2, 1);
			//Matrix2 CC = new Matrix2(new double[,] {{8, 7, 6, 5},{4, 3, 2, 1}});
			//return "AA = " + AA.ToString() + ", BB = " + BB.ToString() + ", CC = " + CC.ToString()
			//	+ "\n, BB*AA = " + (BB.T() * AA.T()).ToString()
			//	+ "\n, DET = " + (BB.T() * AA.T()).det();

			//Swap variables
			//int a=43, b=78;
			//b = (a = a ^ b) ^ b; a = a ^ b;
			//return a.ToString() + ", " + b.ToString();

			//Complex a = 3 + 2 * Complex.i;
			//return a.ToString();

			//i256 a = -1000.0; 
			//i256 a = new i256(false, 1, 2, 3, 4, 5, 6, 7, 8);
			i256 a = 60000, b = 6000;
			return a.ToString() + "\n" + a.ToString(true);

		}

		public static string hex(byte N) { return hex(N, 8); }
		public static string hex(int N) { return hex(N, 32); }
		public static string hex(long N, int m = 64) {
			string s = ""; ulong n; if(N >= 0) { n = (ulong)N; } else { n = (ulong)(~N + 1); }
			for(int q = 0 ; q < m ; q++) { if(n % 2 == 0) { s += "0"; } else { s += "1"; } n >>= 1; }
			return s;
		}
		public static string hex(double N) {
			long n = BitConverter.DoubleToInt64Bits(N);
			int sign = (((n & ~0x7FFFFFFFFFFFFFFFL) >> 63) == 0) ? 1 : -1;
			int exp = (int)((n&0x7FF0000000000000L) >> 52) - 0x3FF;
			long man =     n & 0x000FFFFFFFFFFFFFL | 0x0010000000000000L;
			//                   4030400000000000


			for(int q = 0 ; q < 65 ; q++) { if(man % 2 == 0) { man >>= 1; } else { break; } }

				//6666 5555 5555 5544 4444 4444 3333 3333 3322 2222 2222 1111 1111 1100 0000 0000
				//3210 9876 5432 1098 7654 3210 9876 5432 1098 7654 3210 9876 5432 1098 7654 3210
				//SEEE EEEE EEEE MMMM MMMM MMMM MMMM MMMM MMMM MMMM MMMM MMMM MMMM MMMM MMMM MMMM

				return "sign " + sign.ToString() + ", exp " + exp.ToString() + ", man " + man.ToString("X");
		}

		public static string pfact(int n) {
			List<int> l = listpfact(n); string s = "";
			for(int q = 0 ; q < l.Count - 1 ; q++)
				s += l[q] + ", "; return s + l[l.Count - 1];
		}
		public static List<int> listpfact(int n) {
			List<int> l = new List<int>(); if(n == 0 || n == 1) { l.Add(n); return l; }
			if(n < 0) { n = -n; l.Add(-1); } int i = 2;
			while(n >= 2 && i <= Math.Sqrt(n)) { if(n / (double)i == n / i)
			{ l.Add(i); n /= i; i = 2; } else i++; }
			if(n > 1) l.Add(n); return l;
		}


		public static List<List<int>> perms(List<int> ns) {
			List<List<int>> l = new List<List<int>>(); if(ns.Count == 0 || ns.Count == 1) { l.Add(ns); return l; }


			return l;
		}

		/// <summary>
		/// Discrete Fourier Transform
		/// 
		///https://code.google.com/p/aforge/source/browse/trunk/Sources/Math/FourierTransform.cs
		/// </summary>
		/// <param name="dat">Data array, f()</param>
		/// <param name="f2t">Frequency to Time? or the reverse?</param>
		/// <returns>Transformed data array, F()</returns>
		public static Complex[] DFT(Complex[] dat, bool f2t = true) {
			int n = dat.Length;
			double arg, cos, sin;
			Complex[] dst = new Complex[n];

			for(int i=0 ; i < n ; i++) {
				dst[i] = new Complex(0, 0);



			}

				return dst;
		}
	}

	public class Physics {
		public const double PI = 3.14159265358979323846264338327950288419716939937510; // Diameter to Circumference
		public const double AR = PI / 180.0; // Degrees to Radians
		public const int    c  = 299792458; // Speed of light, in m/s
		public const double G  = 6.673848e-11; // Gravitational Constant, in m3/kg/s2 or N(m/kg)2
		public const double h  = 6.6260695729e-34; // Planck Constant, in J*s or N*m*s or kg*m2/s
		public const double hb = h / 2 / PI; // Reduced Planck Constant, in J*s or N*m*s or kg*m2/s
		public const double kB = 1.380648813e-23; // Bolzmann's Constant, in J/K
		public const double SB = 5.67037321e-8; // Stefan-Boltzmann Constant, in W/m2/K4 or J/s/m2/K4


		/// <summary>
		/// Force of Gravity between two masses at a certain distance
		/// </summary>
		/// <param name="m1">Mass 1</param>
		/// <param name="m2">Mass 2</param>
		/// <param name="r">Distance between them</param>
		/// <returns></returns>
		public static double Fg(double m1, double m2, double r) { return G * ((m1 * m2) / (r * r)); }

	}

	#region Custom Datatypes

	public class i256 : IComparable<i256> {
		public i256() { }
		public i256(bool neg = false, params byte[] ns) { ns.CopyTo(_b, 0); }
		public i256(long n) { BitConverter.GetBytes(n).CopyTo(_b, 0); }
		public i256(double n) { BitConverter.GetBytes(n).CopyTo(_b, 0); }
		public override bool Equals(object obj) { return base.Equals(obj); }
		public override int GetHashCode() { return base.GetHashCode(); }
		public override string ToString() {
			string s = "{";
			for(int q = 0 ; q < 31 ; q++)
				s += _b[q] + (q % 4 == 3 ? "\n" : ", ");
			return s + _b[31] + "}";
		}
		public string ToString(bool b) {
			string s = "{";
			for(int q = 0 ; q < 31 ; q++)
				s += (b ? KN.hex(_b[q]) : _b[q].ToString()) + (q % 4 == 3 ? "\n" : ", ");
			return s + (b ? KN.hex(_b[31]) : _b[31].ToString()) + "}";
		}


		private byte[] _b = new byte[32];// { 0, 0, 0, 0, 0, 0, 0, 0 };


		public static implicit operator long(i256 l) { return BitConverter.ToInt64(l._b, 0); }
		public static implicit operator i256(long l) { return new i256(l); }
		public static implicit operator i256(double l) { return new i256(l); }
		public static bool operator ==(i256 l, i256 r) { for(int q = 0 ; q < 32 ; q++) if(l._b[q] != r._b[q]) return false; return true; }
		public static bool operator !=(i256 l, i256 r) { for(int q = 0 ; q < 32 ; q++) if(l._b[q] != r._b[q]) return true; return false; }
		public static i256 operator ~(i256 l) {
			i256 e = 0;
			for(int q = 0 ; q < 32 ; q++) { e._b[q] = (byte)(~l._b[q] & 255); }
			return e;
		}
		public static i256 operator +(i256 l, i256 r) {
			int ov = 0, v = 0; i256 e = 0;
			for(int q = 0 ; q < 32 ; q++) {
				v = l._b[q] + r._b[q] + ov;
				ov = v >> 8; e._b[q] = (byte)(v & 255);
			}
			return e;
		}



		public int CompareTo(i256 R) { return Sort(this, R); }
		public static int Sort(i256 L, i256 R) {
			if(L.Equals(null) && R.Equals(null)) return 0;
			if(L.Equals(null)) return -1; if(R.Equals(null)) return 1;
			if(L._b[0] == R._b[0]) return 0;
			return 0;
		}

	}
	public class f256 : IComparable<f256> {
		public f256() { }
		//public f256(double n) { _i[0] = (uint)n; }
		public f256(params double[] ns) { }
		public f256(params f256[] ns) {


		}

		private uint[] _i = { 0, 0, 0, 0, 0, 0, 0, 0 };
		//private i256 _i = 0;



		public int CompareTo(f256 R) { return Sort(this, R); }
		public static int Sort(f256 L, f256 R) {
			if(L.Equals(null) && R.Equals(null)) return 0;
			if(L.Equals(null)) return -1; if(R.Equals(null)) return 1;
			if(L._i[0] == R._i[0]) return 0;
			return 0;
		}
	}
	public class fX : IComparable<fX> {
		private uint[] _i = { 0 };
		public int CompareTo(fX R) { return Sort(this, R); }
		public static int Sort(fX L, fX R) {
			if(L.Equals(null) && R.Equals(null)) return 0;
			if(L.Equals(null)) return -1; if(R.Equals(null)) return 1;
			if(L._i[0] == R._i[0]) return 0;
			return 0;
		}
	}

	public class Complex : IComparable<Complex> {
		public Complex(Complex l) { abi(l._a, l._b); }
		public Complex(double na, double nb) { abi(na, nb); }
		public override bool Equals(object obj) { return base.Equals(obj); }
		public override int GetHashCode() { return base.GetHashCode(); }
		public override string ToString() { return ToAbi() + "|" + ToCis(); }
		public string ToAbi() { return "(" + _a.ToString() + (_b >= 0 ? " + " : " - ") + Math.Abs(_b).ToString() + "i)"; }
		public string ToCis() { return "(" + _r.ToString() + ", " + (_t / pi).ToString() + "π)"; }
		public string ToString(string s) { return ToAbi(s) + "|" + ToCis(s); }
		public string ToAbi(string s) { return "(" + _a.ToString(s) + (_b >= 0 ? " + " : " - ") + Math.Abs(_b).ToString(s) + "i)"; }
		public string ToCis(string s) { return "(" + _r.ToString(s) + ", " + (_t / pi).ToString(s) + "π)"; }

		public const double pi = Math.PI;
		public const double tau = 2 * pi;
		public static readonly Complex i = new Complex(0, 1);

		private double _a, _b, _r, _t, _r2;
		public double a { get { return _a; } set { abi(value, _b); } }
		public double b { get { return _b; } set { abi(_a, value); } }
		public double r { get { return _r; } set { cis(value, _t); } }
		public double t { get { return _t; } set { cis(_r, value); } }
		public double r2 { get { return _r2; } set { cis(Math.Sqrt(value), _t); _r2 = value; } }
		public void abi(double na, double nb) {
			_a = na; _b = nb; //_r2 = _a * _a + _b * _b; _r = Math.Sqrt(_r2);
			//if(Math.Abs(_a) < _r / (2 << 24)) _a = 0.0; if(Math.Abs(_b) < _r / (2 << 24)) _b = 0.0;
			_r2 = _a * _a + _b * _b; _r = Math.Sqrt(_r2); _t = Math.Atan2(_b, _a);
		}
		public void cis(double nr, double nth) {
			_r = nr; _t = nth; _a = _r * Math.Cos(_t); _b = _r * Math.Sin(_t); _r2 = _r * _r;
			//if(Math.Abs(_a) < _r / (2 << 24)) _a = 0.0; if(Math.Abs(_b) < _r / (2 << 24)) _b = 0.0;
		}
		public int c {
			get {
				double ca = 1.0, cr = 0.0, cg = 0.0, cb = 0.0;
				double h = (_t + tau) % tau, l = 0.0, d = 0.0;
				//l = ((_r - 1.0) % 1.0 + 1.0) / 1.0;
				//l = ((_r - 0.5) % 0.5 + 0.5) / 1.0;
				//l = ((_r - pi/2.0) % (pi/2.0) + pi/2.0) / pi;
				//l = ((_r - pi) % pi + pi) / tau;
				//d = 0.0;
				if(_r > 16.0) return ~0xFFFFFF + 0xFFFFFF;
				l = _r % 1.0;
				if(Math.Abs(_r) < 1.0 / 128.0) return ~0xFFFFFF; // 0
				if(Math.Abs(_r - Math.Round(_r)) < 1.0 / 128.0) return ~0xFFFFFF + 0xFFFFFF; // Isolines
				//if(h > pi)
				l = Math.Min(_r, 1.0);
				d = Math.Floor(_r)/255.0;
				double H = h * 6.0 / tau; double C = l;// 1.0 - Math.Abs(2.0 * l - 1.0); 
				double X = C * (1.0 - Math.Abs((H % 2.0) - 1.0)), m = l - C;// / 2;
				switch((int)Math.Floor(H)) {
					case 0: cr = C; cg = X; cb = d; break;
					case 1: cr = X; cg = C; cb = d; break;
					case 2: cr = d; cg = C; cb = X; break;
					case 3: cr = d; cg = X; cb = C; break;
					case 4: cr = X; cg = d; cb = C; break;
					case 5: cr = C; cg = d; cb = X; break;
				} return ~0xFFFFFF +
					((int)Math.Floor((cr + m) * 255.0) << 16) +
					((int)Math.Floor((cg + m) * 255.0) << 8) +
					 (int)Math.Floor((cb + m) * 255.0);
			}
			set {
				int ca, cr, cg, cb; byte[] cc = new byte[4];
				cc = BitConverter.GetBytes(value);
				ca = cc[0]; cr = cc[1]; cg = cc[2]; cb = cc[3];
				double al = (2.0 * cr - cg - cb) / 2.0, be = (cg - cb) / 2.0 * Math.Sqrt(3);
				abi(al, be);



			}
		}

		public Complex conj() { return new Complex(_a, -_b); }
		public Complex norm() { return this / _r; }
		public Complex sqrt() { return new Complex(Math.Sqrt((_a + _r) / 2), Math.Sign(_b) * Math.Sqrt((-_a + _r) / 2)); }
		public Complex ln() { return new Complex(Math.Log(_r), _t); }
		public Complex exp() { return Cis(Math.Exp(_a), _b); }

		public Complex sin() { return new Complex(Math.Sin(_a) * Math.Cosh(_b), Math.Cos(_a) * Math.Sinh(_b)); }
		public Complex cos() { return new Complex(Math.Cos(_a) * Math.Cosh(_b), -Math.Sin(_a) * Math.Sinh(_b)); }
		public Complex tan() { return sin() / cos(); }

		public static implicit operator Complex(double l) { return new Complex(l, 0); }
		public static bool operator ==(Complex l, Complex r) { return (l._a == r._a) && (l._b == r._b); }
		public static bool operator !=(Complex l, Complex r) { return !((l._a == r._a) && (l._b == r._b)); }
		public static double dot(Complex l, Complex r) { return l._a * r._a + l._b * r._b; }
		public static Complex Cis(double r, double th) { return new Complex(r * Math.Cos(th), r * Math.Sin(th)); }
		public static Complex operator +(Complex l, Complex r) { return new Complex(l._a + r._a, l._b + r._b); }
		public static Complex operator +(Complex l, double r) { return new Complex(l._a + r, l._b); }
		public static Complex operator +(double l, Complex r) { return new Complex(l + r._a, r._b); }
		public static Complex operator -(Complex l, Complex r) { return new Complex(l._a - r._a, l._b - r._b); }
		public static Complex operator -(Complex l, double r) { return new Complex(l._a - r, l._b); }
		public static Complex operator -(double l, Complex r) { return new Complex(l - r._a, -r._b); }
		public static Complex operator *(Complex l, Complex r) {
			return new Complex(l._a * r._a - l._b * r._b, l._a * r._b + r._a * l._b);
		}
		public static Complex operator *(Complex l, double r) { return new Complex(l._a * r, l._b * r); }
		public static Complex operator *(double l, Complex r) { return new Complex(l * r._a, l * r._b); }
		public static Complex operator /(Complex l, Complex r) {
			return r._r2 == 0 ? new Complex(double.NaN, double.NaN) :
			new Complex((l._a * r._a + l._b * r._b) / r._r2, (l._b * r._a - l._a * r._b) / r._r2);
		}
		public static Complex operator /(Complex l, double r) {
			return r == 0 ? new Complex(double.NaN, double.NaN) : new Complex(l._a / r, l._b / r);
		}
		public static Complex operator /(double l, Complex r) {
			return r._r2 == 0 ? new Complex(double.NaN, double.NaN) :
			new Complex((l * r._a) / r._r2, (-l * r._b) / r._r2);
		}
		public static Complex operator ^(Complex l, Complex r) {
			if(r._b == 0) return l ^ r._a;
			return Cis(Math.Pow(l._r, r._a) * Math.Exp(-r._b * l._t), r._b * Math.Log(l._r) + r._a * l._t);
		}
		public static Complex operator ^(Complex l, double r) { return Cis(Math.Pow(l._r, r), r * l._t); }
		public static Complex operator ^(double l, Complex r) {
			double lr = Math.Abs(l), lt = l < 0 ? Math.PI : 0;
			return Cis(Math.Pow(lr, r._a) * Math.Exp(-r._b * lt), r._b * Math.Log(Math.Abs(l)) + r._a * lt);
		}


		public int CompareTo(Complex R) { return SortR(this, R); }
		public static int SortA(Complex L, Complex R) {
			if(L.Equals(null) && R.Equals(null)) return 0;
			if(L.Equals(null)) return -1; if(R.Equals(null)) return 1;
			if(L._a == R._a) return L._b.CompareTo(R._b); return L._a.CompareTo(R._a);
		}
		public static int SortB(Complex L, Complex R) {
			if(L.Equals(null) && R.Equals(null)) return 0;
			if(L.Equals(null)) return -1; if(R.Equals(null)) return 1;
			if(L._b == R._b) return L._a.CompareTo(R._a); return L._b.CompareTo(R._b);
		}
		public static int SortR(Complex L, Complex R) {
			if(L.Equals(null) && R.Equals(null)) return 0;
			if(L.Equals(null)) return -1; if(R.Equals(null)) return 1;
			if(L._r == R._r) return L._t.CompareTo(R._t); return L._r.CompareTo(R._r);
		}
		public static int SortT(Complex L, Complex R) {
			if(L.Equals(null) && R.Equals(null)) return 0;
			if(L.Equals(null)) return -1; if(R.Equals(null)) return 1;
			if(L._t == R._t) return L._r.CompareTo(R._r); return L._t.CompareTo(R._t);
		}


	}

	public class Quaternion : IComparable<Quaternion> {
		public Quaternion(Quaternion l) { abi(l._a, l._b, l._c, l._d); }
		public Quaternion(Complex na, Complex nb) { abi(na.a, na.b, nb.a, nb.b); }
		public Quaternion(double na, double nb, double nc, double nd) { abi(na, nb, nc, nd); }
		public override bool Equals(object obj) { return base.Equals(obj); }
		public override int GetHashCode() { return base.GetHashCode(); }
		public override string ToString() { return ToAbi() + "|" + ToCis(); }
		public string ToAbi() {
			return "(" + _a.ToString() +
				(_b >= 0 ? " + " : " - ") + Math.Abs(_b).ToString() + "i" +
				(_c >= 0 ? " + " : " - ") + Math.Abs(_c).ToString() + "j" +
				(_d >= 0 ? " + " : " - ") + Math.Abs(_d).ToString() + "k)";
		}
		public string ToCis() { return "(" + _r.ToString() + ", " + V.norm().ToAbi() + ", " + _t.ToString() + "π)"; }

		public const double pi = Math.PI;
		public const double tau = 2 * pi;
		public static readonly Quaternion i = new Quaternion(0, 1, 0, 0);
		public static readonly Quaternion j = new Quaternion(0, 0, 1, 0);
		public static readonly Quaternion k = new Quaternion(0, 0, 0, 1);

		private double _a, _b, _c, _d, _r, _t, _r2, _v, _v2;
		public double a { get { return _a; } set { abi(value, _b, _c, _d); } }
		public double b { get { return _b; } set { abi(_a, value, _c, _d); } }
		public double c { get { return _c; } set { abi(_a, _b, value, _d); } }
		public double d { get { return _d; } set { abi(_a, _b, _c, value); } }
		public double r { get { return _r; } set { cis(value, V, _t); } }
		public double t { get { return _t; } set { cis(_r, V, value); } }
		public double r2 { get { return _r2; } set { cis(Math.Sqrt(value), V, _t); _r2 = value; } }
		public double v { get { return _v; } }
		public double v2 { get { return _v2; } }
		public Quaternion V {
			get { return new Quaternion(0, _b, _c, _d); }
			set { value = value.norm(); abi(_a, value._b, value._c, value._d); }
		}
		public void abi(double na, double nb, double nc, double nd) {
			_a = na; _b = nb; _c = nc; _d = nd; _v2 = _b * _b + _c * _c + _d * _d;
			_v = Math.Sqrt(_v2); _r2 = _a * _a + _v2; _r = Math.Sqrt(_r2); _t = Math.Acos(_a / _r);
		}
		public void cis(double nr, Quaternion nv, double nth) {
			nv = nv.norm();
			_r = nr;
			_t = nth;
			_a = _r * Math.Cos(_t);
			_b = _r * nv._b * Math.Sin(_t);
			_c = _r * nv._c * Math.Sin(_t);
			_d = _r * nv._d * Math.Sin(_t);
			_r2 = _r * _r;
			_v2 = _r2 - _a * _a;
			_v = Math.Sqrt(_v2);
		}


		public Quaternion conj() { return new Quaternion(_a, -_b, -_c, -_d); }
		public Quaternion norm() { return this / _r; }
		//public Complex sqrt() { return new Complex(Math.Sqrt((_a + _r) / 2), Math.Sign(_b) * Math.Sqrt((-_a + _r) / 2)); }
		public Quaternion ln() { return Math.Log(_r) + V.norm() * Math.Acos(_a / _r); }
		public Quaternion exp() { return Math.Exp(_a) * (Math.Cos(_v) + V.norm() * Math.Sin(_v)); }

		public static bool operator ==(Quaternion l, Quaternion r) {
			return (l._a == r._a) && (l._b == r._b) && (l._c == r._c) && (l._d == r._d);
		}
		public static bool operator !=(Quaternion l, Quaternion r) {
			return !((l._a == r._a) && (l._b == r._b) && (l._c == r._c) && (l._d == r._d));
		}
		public static double dot(Quaternion l, Quaternion r) { return l._a * r._a + l._b * r._b + l._c * r._c + l._d * r._d; }
		public static Quaternion Cis(double r, Quaternion n, double th) {
			n = r * n.norm() * Math.Sin(th);
			return new Quaternion(r * Math.Cos(th), n.b, n.c, n.d);
		}
		public static Quaternion Cis(double r, double b, double c, double d, double th) {
			return Cis(r, new Quaternion(0, b, c, d), th);
		}
		public static Quaternion operator +(Quaternion l, Quaternion r) {
			return new Quaternion(l._a + r._a, l._b + r._b, l._c + r._c, l._d + r._d);
		}
		public static Quaternion operator +(Quaternion l, double r) { return new Quaternion(l._a + r, l._b, l._c, l._d); }
		public static Quaternion operator +(double l, Quaternion r) { return new Quaternion(l + r._a, r._b, r._c, r._d); }
		public static Quaternion operator -(Quaternion l, Quaternion r) {
			return new Quaternion(l._a - r._a, l._b - r._b, l._c - r._c, l._d - r._d);
		}
		public static Quaternion operator -(Quaternion l, double r) { return new Quaternion(l._a - r, l._b, l._c, l._d); }
		public static Quaternion operator -(double l, Quaternion r) { return new Quaternion(l - r._a, -r._b, -r._c, -r._d); }
		public static Quaternion operator *(Quaternion l, Quaternion r) {
			return new Quaternion(
				l._a * r._a - l._b * r._b - l._c * r._c - l._d * r._d, l._a * r._b + l._b * r._a + l._c * r._d - l._d * r._c,
				l._a * r._c - l._b * r._d + l._c * r._a + l._d * r._b, l._a * r._d + l._b * r._c - l._c * r._b + l._d * r._a);
		}
		public static Quaternion operator *(Quaternion l, double r) { return new Quaternion(l._a * r, l._b * r, l._c * r, l._d * r); }
		public static Quaternion operator *(double l, Quaternion r) { return new Quaternion(l * r._a, l * r._b, l * r._c, l * r._d); }
		public static Quaternion operator /(Quaternion l, Quaternion r) {
			return r._r2 == 0 ?
				new Quaternion(double.NaN, double.NaN, double.NaN, double.NaN) : l * (1 / r);
		}
		public static Quaternion operator /(Quaternion l, double r) {
			return r == 0 ? new Quaternion(double.NaN, double.NaN, double.NaN, double.NaN)
				: new Quaternion(l._a / r, l._b / r, l._c / r, l._d / r);
		}
		public static Quaternion operator /(double l, Quaternion r) {
			return r._r2 == 0 ? new Quaternion(double.NaN, double.NaN, double.NaN, double.NaN) :
			new Quaternion((l * r._a) / r._r2, (-l * r._b) / r._r2, (-l * r._c) / r._r2, (-l * r._d) / r._r2);
		}
		//public static Complex operator ^(Complex l, Complex r) {
		//	if(r._b == 0) return l ^ r._a;
		//	return Cis(Math.Pow(l._r, r._a) * Math.Exp(-r._b * l._t), r._b * Math.Log(l._r) + r._a * l._t);
		//}
		//public static Complex operator ^(Complex l, double r) { return Cis(Math.Pow(l._r, r), r * l._t); }
		//public static Complex operator ^(double l, Complex r) {
		//	double lr = Math.Abs(l), lt = l < 0 ? Math.PI : 0;
		//	return Cis(Math.Pow(lr, r._a) * Math.Exp(-r._b * lt), r._b * Math.Log(Math.Abs(l)) + r._a * lt);
		//}


		public int CompareTo(Quaternion R) { return SortR(this, R); }
		public static int SortA(Quaternion L, Quaternion R) {
			if(L.Equals(null) && R.Equals(null)) return 0;
			if(L.Equals(null)) return -1; if(R.Equals(null)) return 1;
			if(L._a == R._a) return L._b.CompareTo(R._b); return L._a.CompareTo(R._a);
		}
		public static int SortB(Quaternion L, Quaternion R) {
			if(L.Equals(null) && R.Equals(null)) return 0;
			if(L.Equals(null)) return -1; if(R.Equals(null)) return 1;
			if(L._b == R._b) return L._a.CompareTo(R._a); return L._b.CompareTo(R._b);
		}
		public static int SortR(Quaternion L, Quaternion R) {
			if(L.Equals(null) && R.Equals(null)) return 0;
			if(L.Equals(null)) return -1; if(R.Equals(null)) return 1;
			if(L._r == R._r) return L._t.CompareTo(R._t); return L._r.CompareTo(R._r);
		}
		public static int SortT(Quaternion L, Quaternion R) {
			if(L.Equals(null) && R.Equals(null)) return 0;
			if(L.Equals(null)) return -1; if(R.Equals(null)) return 1;
			if(L._t == R._t) return L._r.CompareTo(R._r); return L._t.CompareTo(R._t);
		}


	}

	public class Matrix2 {
		public Matrix2() : this(1, 1, 0.0) { }
		public Matrix2(int nx) {
			_x = nx; _y = nx;
			_dat = new double[nx,nx]; for(int y = 0 ; y < nx ; y++)
				for(int x = 0 ; x < nx ; x++) { _dat[y,x] = (x == y) ? 1 : 0; }
		}
		public Matrix2(Matrix2 l) : this(l._dat) { }
		public Matrix2(int nx, int ny, params double[] nd) {
			_x = nx; _y = ny; int le = nd.GetLength(0), lo = 0;
			_dat = new double[ny,nx]; for(int y = 0 ; y < ny ; y++)
				for(int x = 0 ; x < nx ; x++) { lo = y * nx + x; _dat[y,x] = (lo < le) ? nd[lo] : 0; }
		}
		public Matrix2(double[,] nd) {
			_y = nd.GetLength(0); _x = nd.GetLength(1);
			_dat = new double[_y,_x]; for(int y = 0 ; y < _y ; y++)
				for(int x = 0 ; x < _x ; x++) { _dat[y,x] = nd[y,x]; }
		}
		public Matrix2(double[] nd, bool row = false) {
			if(row) { _y = 1; _x = nd.GetLength(0);
				_dat = new double[_y, _x]; for(int x = 0 ; x < _x ; x++) { _dat[0, x] = nd[x]; }
			} else { _y = nd.GetLength(0); _x = 1;
				_dat = new double[_y, _x]; for(int y = 0 ; y < _y ; y++) { _dat[y, 0] = nd[y]; }
			}
		}

		public override bool Equals(object obj) { return base.Equals(obj); }
		public override int GetHashCode() { return base.GetHashCode(); }
		public override string ToString() {
			string s = "{"; for(int y = 0 ; y < _dat.GetLength(0) ; y++) { s += "[ ";
				for(int x = 0 ; x < _dat.GetLength(1) ; x++) { s += _dat[y,x] + " "; } s += "]";
			} return (s += "}");
		}

		private double[,] _dat; private readonly int _x, _y;
		public double this[int x, int y] { get { return _dat[y,x]; } set { _dat[y,x] = value; } }

		public Matrix2 T() {
			Matrix2 t = new Matrix2(_y, _x); for(int y = 0 ; y < _y ; y++) {
				for(int x = 0 ; x < _x ; x++) { t._dat[x,y] = _dat[y,x]; }
			} return t;
		}
		public Matrix2 row(int y) { Matrix2 t = new Matrix2(_x, 1);
			for(int x = 0 ; x < _x ; x++) { t._dat[0, x] = _dat[y, x]; } return t; }
		public Matrix2 col(int x) { Matrix2 t = new Matrix2(1, _y);
			for(int y = 0 ; y < _y ; y++) { t._dat[y, 0] = _dat[y, x]; } return t; }
		
		public Matrix2 sub(int x1, int y1, int x2, int y2) {
			Matrix2 t = new Matrix2(x2 - x1, y2 - y1);
			for(int y = y1 ; y < y2 ; y++) {
				for(int x = x1 ; x < x2 ; x++) {
					if(x > _x || y > _y) t._dat[y,x] = 0.0;
					else t._dat[y,x] = _dat[y - y1,x - x1];
				}
			} return t;
		}
		public Matrix2 minor(int x1, int y1) {
			//if(x1 < 0 || x1 >= _x || y1 < 0 || y1 >= _y) return this;
			x1 = Math.Max(Math.Min(x1, _x - 1), 0); y1 = Math.Max(Math.Min(y1, _y - 1), 0);
			Matrix2 t = new Matrix2(_x - 1, _y - 1);
			for(int y = 0 ; y < _y - 1 ; y++) {//if(y == y1) continue;
				for(int x = 0 ; x < _x - 1 ; x++) {//if(x == x1) continue;
					t._dat[y,x] = _dat[y + (y >= y1 ? 1 : 0),x + (x >= x1 ? 1 : 0)];
				}
			}
			return t;
		}
		public double det() {
			if(_x != _y) return Double.NaN;
			if(_x == 1 && _y == 1) return _dat[0,0];
			double d = 0, dr = 0, dl = 0;
			Matrix2 m;
			for(int x = 0 ; x < _x ; x++) {
				m = minor(x, 0); dr = _dat[0,x]; dl = (x % 2 == 1 ? -1 : 1);
				d += dl * dr * m.det();

			}
			return d;
		}

		public static Matrix2 operator +(Matrix2 l, Matrix2 r) {
			Matrix2 t = new Matrix2(Math.Max(l._y, r._y), Math.Max(l._x, r._x));
			for(int y = 0 ; y < t._y ; y++) {
				for(int x = 0 ; x < t._x ; x++) {
					t._dat[y,x] = ((y < l._y && x < l._x) ? l._dat[y,x] : 0)
						+ ((y < r._y && x < r._x) ? r._dat[y,x] : 0);
				}
			} return t;
		}
		public static Matrix2 operator +(Matrix2 l, double r) {
			Matrix2 t = new Matrix2(l._y, l._x);
			for(int y = 0 ; y < t._y ; y++) {
				for(int x = 0 ; x < t._x ; x++) {
					t._dat[y,x] = l._dat[y,x] + r;
				}
			} return t;
		}
		public static Matrix2 operator +(double l, Matrix2 r) {
			Matrix2 t = new Matrix2(r._y, r._x);
			for(int y = 0 ; y < t._y ; y++) {
				for(int x = 0 ; x < t._x ; x++) {
					t._dat[y,x] = l + r._dat[y,x];
				}
			} return t;
		}
		public static Matrix2 operator -(Matrix2 l, Matrix2 r) {
			Matrix2 t = new Matrix2(Math.Max(l._y, r._y), Math.Max(l._x, r._x));
			for(int y = 0 ; y < t._y ; y++) {
				for(int x = 0 ; x < t._x ; x++) {
					t._dat[y,x] = ((y < l._y && x < l._x) ? l._dat[y,x] : 0)
						- ((y < r._y && x < r._x) ? r._dat[y,x] : 0);
				}
			} return t;
		}
		public static Matrix2 operator -(Matrix2 l, double r) {
			Matrix2 t = new Matrix2(l._y, l._x);
			for(int y = 0 ; y < t._y ; y++) {
				for(int x = 0 ; x < t._x ; x++) {
					t._dat[y,x] = l._dat[y,x] - r;
				}
			} return t;
		}
		public static Matrix2 operator -(double l, Matrix2 r) {
			Matrix2 t = new Matrix2(r._y, r._x);
			for(int y = 0 ; y < t._y ; y++) {
				for(int x = 0 ; x < t._x ; x++) {
					t._dat[y,x] = l - r._dat[y,x];
				}
			} return t;
		}
		//public static Matrix2 operator *(Matrix2 l, Matrix2 r) {
		//	Matrix2 t = new Matrix2(Math.Max(l._y, r._y), Math.Max(l._x, r._x));
		//	int z = 0;
		//	for(int y = 0 ; y < t._y ; y++) {
		//		for(int x = 0 ; x < t._x ; x++) {

		//			t._dat[y,x] = ((y < l._y && x < l._x) ? l._dat[y,x] : 0)
		//				- ((y < r._y && x < r._x) ? r._dat[y,x] : 0);
		//		}
		//	} return t;
		//}
		public static Matrix2 operator *(Matrix2 l, Matrix2 r) {
			if(l._y != r._x) return new Matrix2(1, 1, Double.NaN);
			Matrix2 t = new Matrix2(r._x, l._y);
			for(int y = 0 ; y < t._y ; y++) {
				for(int x = 0 ; x < t._x ; x++) {
					for(int z = 0 ; z < l._x ; z++) {
						t._dat[y, x] += l._dat[y, z] * r._dat[z, x];
					}
				}
			} return t;
		}
		public static Matrix2 operator *(Matrix2 l, double r) {
			Matrix2 t = new Matrix2(l._y, l._x);
			for(int y = 0 ; y < t._y ; y++) {
				for(int x = 0 ; x < t._x ; x++) {
					t._dat[y,x] = l._dat[y,x] * r;
				}
			} return t;
		}
		public static Matrix2 operator *(double l, Matrix2 r) {
			Matrix2 t = new Matrix2(r._y, r._x);
			for(int y = 0 ; y < t._y ; y++) {
				for(int x = 0 ; x < t._x ; x++) {
					t._dat[y,x] = l * r._dat[y,x];
				}
			} return t;
		}




	}

	#endregion Custom Datatypes


	public class NeuralNetwork {
		public struct Neuron {
			public int numinputs; public List<double> weights;
			public Neuron(int nin) {
				numinputs = nin + 1; weights = new List<double>(); Random rnd = new Random();
				for(int q = 0 ; q < numinputs ; q++) { weights.Add(rnd.NextDouble()); }
			}
		}
		public struct Layer {
			public int numneurons; public List<Neuron> neurons;
			public Layer(int nin, int ninn) {
				numneurons = nin; neurons = new List<Neuron>();
				for(int q = 0 ; q < numneurons ; q++) { neurons.Add(new Neuron(ninn)); }
			}
		}

		private int numinputs, numoutputs, numlayers, numperlayer;
		private List<Layer> layers;

		public NeuralNetwork() {


		}
		public void CreateNet() {


		}
		public List<double> GetWeights() {
			List<double> w = new List<double>();
			for(int x = 0 ; x < numlayers ; x++ ) for(int y = 0 ; y < layers[x].numneurons ; y++)
			for(int z = 0 ; z < layers[x].neurons[y].numinputs ; z++) w.Add(layers[x].neurons[y].weights[z]);
			return w;
		}
		public void PutWeights(List<double> w) { int p = 0;
			for(int x = 0 ; x < numlayers ; x++) for(int y = 0 ; y < layers[x].numneurons ; y++)
			for(int z = 0 ; z < layers[x].neurons[y].numinputs ; z++) layers[x].neurons[y].weights[z] = w[p++];
		}

		public List<double> Update(List<double> inputs) {
			List<double> outputs = new List<double>(); int cweight = 0;
			if(inputs.Count != numinputs) return outputs;
			double netinput = 0; int numinp = 0;

			for(int x = 0 ; x < numlayers ; x++) {
				if(x > 0) inputs = outputs; outputs.Clear();
				for(int y = 0 ; y < layers[x].numneurons ; y++) {
					cweight = 0; netinput = 0; numinp = layers[x].neurons[y].numinputs;
					for(int z = 0 ; z < numinp ; z++) {
						netinput += layers[x].neurons[y].weights[z] * ((z == 0) ? -1.0 : inputs[cweight++]);
					}
					outputs.Add(sig(netinput));
				}
			}
			return outputs;
		}


		public static double sig(double a, double p = 1.0) { return 1.0 / (1.0 + Math.Exp(-a / p)); }
	}
	public class GenAlg {




	}


}
