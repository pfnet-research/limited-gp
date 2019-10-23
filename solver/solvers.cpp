#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <queue>
#include <fstream>
#include <sstream>
using namespace std;
class XorShift {
private:
	
	unsigned int x, y, z, w;

public:

	XorShift() {
		init();
	}

	void init() {
		x = 314159265; y = 358979323; z = 846264338; w = 327950288;
	}

	void setSeed(int seed) {
		x = 314159265; y = 358979323; z = 846264338; w = 327950288 + seed;
		for (int i = 0; i < 500; i++) {
			next();
		}
	}

	inline unsigned int next() {
		unsigned int t = x ^ x << 11; x = y; y = z; z = w; return w = w ^ w >> 19 ^ t ^ t >> 8;
	}

	inline int nextInt(int m) {
		return (int)(next() % m);
	}


	inline int nextInt(int min, int max) {
		return min + nextInt(max - min + 1);
	}


	inline double nextDouble() {
		return (double)next() / ((long long)1 << 32);
	}

	inline double nextGaussian() {
		return sqrt(-2 * log(nextDouble())) * cos(2 * 3.14159265358979323846 * nextDouble());
	}

};
XorShift rnd;
int D, N;
double N0I;
int cnt = 0;
vector<double> low;
vector<double> high;
double best_so_far = 1e9;

double ask(const vector<double> &x) {
	cnt++;
	if (cnt > N) {
		cerr << "eval limit exceeded" << endl;
	}

	for (int i = 0; i < D; i++) {
		cout << max(low[i], min(high[i], x[i])) << (i == D - 1 ? '\n' : ' ');
	}
	cout.flush();
	double res;
	cin >> res;
	for (int i = 0; i < D; i++) {
		res += max((double)0.0, 1e5*(x[i] - high[i]));
		res += max((double)0.0, 1e5*(low[i] - x[i]));
	}
	return res;
}

void input() {
	cout << fixed << setprecision(15);
	cerr << fixed << setprecision(15);
	N0I = sqrt(D) * (1 - (1.0 / (4 * D)) + (1.0 / (21 * D * D)));
	cin >> D >> N;
	cnt = 0;
	low.clear(); low.resize(D);
	high.clear(); high.resize(D);
	for (int i = 0; i < D; i++) {
		cin >> low[i];
	}
	for (int i = 0; i < D; i++) {
		cin >> high[i];
	}
	best_so_far = 1e9;
}

vector<double> get_random_point() {
	vector<double> X(D);
	for (int i = 0; i < D; i++) {
		X[i] = (high[i] - low[i]) * rnd.nextDouble() + low[i];
	}
	return X;
}

vector<vector<double> > inverse(const vector<vector<double> > &A) {

	vector<vector<double> > B(A.size(), vector<double>(2 * A.size(), 0));
	for (int i = 0; i < A.size(); i++) {
		for (int j = 0; j < A.size(); j++) {
			B[i][j] = A[i][j];
		}
	}

	for (int i = 0; i < A.size(); i++) {

		B[i][A.size() + i] = 1;

	}

	for (int i = 0; i < A.size(); i++) {
		int pivot = i;
		for (int j = i; j < A.size(); j++) {
			if (abs(B[j][i]) > abs(B[pivot][i])) {
				pivot = j;
			}
		}
		swap(B[i], B[pivot]);
		for (int j = i + 1; j < 2 * A.size(); j++) {
			B[i][j] /= B[i][i];
		}

		for (int j = 0; j < A.size(); j++) {
			if (i != j) {
				for (int k = i + 1; k < 2 * A.size(); k++) {
					B[j][k] -= B[j][i] * B[i][k];
				}
			}
		}
	}
	vector<vector<double> > res(A.size(), vector<double>(A.size()));

	for (int i = 0; i < A.size(); i++) {
		for (int j = 0; j < A.size(); j++) {
			res[i][j] = B[i][A.size() + j];
		}
	}

	return res;

}

vector<vector<double> > cholesky(const vector<vector<double> > &A) {
	vector<vector<double> > L(D, vector<double>(D, 0));

	for (int i = 0; i < D; i++)
		for (int j = 0; j < (i + 1); j++) {
			double s = 0;
			for (int k = 0; k < j; k++)
				s += L[i][k] * L[j][k];
			L[i][j] = (i == j) ?
				sqrt(max(1e-20, A[i][i] - s)) :
				(1.0 / L[j][j] * (A[i][j] - s));
		}

	return L;
}

vector<double> Gaussian(const vector<double> &m, const vector<vector<double> > &sigma, double rate = 1) {
	vector<double> g(D), res(D);
	for (int i = 0; i < D; i++) {
		g[i] = rnd.nextGaussian();
	}
	vector<vector<double> > c_sigma = cholesky(sigma);
	for (int i = 0; i < D; i++) {
		res[i] = m[i];
		for (int j = 0; j < D; j++) {
			res[i] += rate * g[j] * c_sigma[i][j];
		}
	}
	return res;
}

void calc_sigma(const vector<vector<double> >&v,
	vector<double> weight,
	vector<double> &mean,
	vector<vector<double> > &sigma,
	bool fixed_mean = false) {
	sigma.clear(); sigma.resize(D, vector<double>(D, 0));
	if (!fixed_mean) {
		double sum = 0;
		mean.clear(); mean.resize(D);
		for (int i = 0; i < weight.size(); i++)sum += weight[i];
		for (int i = 0; i < weight.size(); i++)weight[i] /= sum;
		for (int i = 0; i < D; i++) {
			mean[i] = 0;
			for (int j = 0; j < v.size(); j++) {
				mean[i] += v[j][i] * weight[j];
			}
		}
	}
	for (int i = 0; i < D; i++) {
		for (int j = 0; j < D; j++) {
			sigma[i][j] = 0;
			for (int k = 0; k < v.size(); k++) {
				sigma[i][j] += weight[k] * (v[k][i] - mean[i]) * (v[k][j] - mean[j]);
			}
		}
	}
}

inline double dst(const vector<double> &A, const vector<double> &B) {
	double res = 0;
	for (int i = 0; i < A.size(); i++) {
		res += (A[i] - B[i]) * (A[i] - B[i]);
	}
	return res;
}

inline double dst_from_origin(const vector<double> &A) {
	double res = 0;
	for (int i = 0; i < A.size(); i++) {
		res += A[i] * A[i];
	}
	return res;
}

void init_mean_sigma(vector<double> &mean, vector<vector<double> > &sigma) {
	mean.clear(); mean.resize(D);
	sigma.clear(); sigma.resize(D, vector<double>(D, 0));
	for (int i = 0; i < D; i++) {
		mean[i] = (low[i] + high[i]) * 0.5;
		sigma[i][i] = 0.04 * (high[i] - low[i]) * (high[i] - low[i]);
	}

}

struct Point {
	vector<double> p;
	double score;
	int idx;
	bool evaled = false;

	Point() {}

	Point(const vector<double> &_p, double _score) {
		p = _p;
		score = _score;
	}

	bool operator<(const Point &right) const {
		return score < right.score;
	}
};

class Limited_GP {
private:

	int NEvals = 0;
	int cur = 0;
	vector<pair<double, vector<double> > > pq;
	vector<double> ys;
	vector<vector<double> > invK;
	vector<double> mean;
	vector<vector<double> > sigma;
	vector<vector<double> > _ss;
	vector<vector<double> > ss;
	vector<double> M;
	vector<Point> vp;
	vector<double> bestS;
	vector<double> curS;
	double ac_coef = 2.0;
	double c1 = 0.95;
	double c2 = 0.05;
	double c4 = 0.1;
	int Z = 15;

public:

	Limited_GP() {
		bestS.resize(1000 * D, 1e18);
		curS = bestS;
		M.resize(D, 0);
		init_mean_sigma(mean, sigma);
		Z = 7 + 3 * log(D);
		c2 = min(0.99, 3.6 / (D * D));
		c1 = 1 - c2;
		c4 = max(0.0, 1 - (3.6 / (D * D)));
	}

	int get_NEvals() {
		return NEvals;
	}

	void update_bestS() {
		bestS = curS;
	}

	pair<double, vector<double> > get_best_point() {
		return pq[0];
	}

	double kernel(const vector<double> &X1, const vector<double> &X2) {
		return exp(-dst(X1, X2));
	}

	void init(vector<double> &_mean, vector<vector<double> > &_sigma) {
		NEvals = 0;
		pq.clear();
		mean = _mean;
		sigma = _sigma;
		M.clear();
		M.resize(D, 0);
	}
	bool converge() {
		if (NEvals > 100 * D) {
			if (bestS[NEvals / 3] < pq[0].first) {
				return true;
			}
		}
		if (NEvals > 100 * D) {
			if (curS[NEvals - 100 * D] - 1e-5 < pq[0].first) {
				return true;
			}
		}
		return (int)pq.size() > D && (pq[D].first - pq[0].first < 1e-9);
	}
	double th;

	void preprocess() {
		vector<vector<double> > v;
		int mx = 1000 * D;
		double prog = (double)NEvals / mx;

		int sz = min((int)pq.size(), (int)(Z * (1 - 0.8 * prog)));
		for (int i = 0; i < sz; i++) {
			v.push_back(pq[i].second);
		}

		vector<double> W(v.size(), 1);
		double r = pow(1 - prog, 1.0 / v.size());
		for (int i = 0; i < v.size(); i++) {
			W[i] = pow(r, i);
		}

		vector<double> nmean;
		vector<vector<double> > nsigma;
		calc_sigma(v, W, nmean, nsigma, false);

		for (int i = 0; i < D; i++) {
			M[i] = c4 * M[i] + 1.0 * (nmean[i] - mean[i]);
			mean[i] = mean[i] + M[i];
		}

		for (int i = 0; i < D; i++) {
			for (int j = 0; j < D; j++) {
				sigma[i][j] = c1 * sigma[i][j] + c2 * nsigma[i][j];
			}
		}

		_ss = cholesky(sigma);
		ss = inverse(_ss);
		vp.clear();
		vp.resize(pq.size());

		for (int i = 0; i < pq.size(); i++) {
			vp[i].p.resize(D);
			vp[i].idx = i;

			for (int j = 0; j < D; j++) {
				vp[i].p[j] = 0;
				for (int k = 0; k < D; k++) {
					vp[i].p[j] += ss[j][k] * (pq[i].second[k] - mean[k]);
				}
			}

			vp[i].score = dst_from_origin(vp[i].p);
			if (vp[i].score > 1e18)vp[i].score = 1e18;
		}

		sort(vp.begin(), vp.end());

		if ((int)vp.size() > 2 * sz)vp.resize(2 * sz);
		th = vp[(int)vp.size() - 1].score;
		th = vp[sz - 1].score;
		invK.clear();
		invK.resize(vp.size(), vector<double>(vp.size(), 0));
		vector<vector<double> > K = invK;

		for (int i = 0; i < vp.size(); i++) {
			for (int j = i; j < vp.size(); j++) {
				K[i][j] = kernel(vp[i].p, vp[j].p);
				K[j][i] = K[i][j];
			}
		}
		invK = inverse(K);
	}

	double calc(const vector<double> &X, double a) {

		vector<double> k(vp.size());
		for (int i = 0; i < k.size(); i++) {
			k[i] = kernel(X, vp[i].p);
		}

		double m = 0;
		for (int i = 0; i < vp.size(); i++) {
			for (int j = 0; j < vp.size(); j++) {
				m += invK[i][j] * vp[j].idx * k[i];
			}
		}

		double s = 0;
		for (int i = 0; i < vp.size(); i++) {
			for (int j = 0; j < vp.size(); j++) {
				s += invK[i][j] * k[j] * k[i];
			}
		}

		s = kernel(X, X) - s;
		return m - a * s;
	}

	vector<double> grad(const vector<double> &X, double a) {
		vector<double> g(D, 0);
		vector<double> k(vp.size());
		for (int i = 0; i < k.size(); i++) {
			k[i] = kernel(X, vp[i].p);
		}
		double s = 0;
		for (int i = 0; i < vp.size(); i++) {
			for (int j = 0; j < vp.size(); j++) {
				s = 2 * invK[i][j] * k[j] * k[i];
				for (int d = 0; d < D; d++) {
					g[d] -= s * ((X[d] - vp[i].p[d]) + (X[d] - vp[j].p[d]));
				}
			}
		}
		for (int i = 0; i < D; i++) {
			g[i] *= a;
		}
		double m = 0;
		for (int i = 0; i < vp.size(); i++) {
			double t = 0;
			for (int j = 0; j < vp.size(); j++) {
				t += invK[i][j] * vp[j].idx;
			}
			t *= k[i];
			for (int j = 0; j < D; j++) {
				g[j] -= 2 * t * (X[j] - vp[i].p[j]);
			}
		}
		return g;
	}

	vector<double> get_next_sample() {
		NEvals++;
		vector<double> X;
		if ((int)pq.size() >= Z) {
			preprocess();
			vector<double> tX(D);
			double bestscore = 1e9;
			vector<double> best;
			for (int loop = 0; loop < 3; loop++) {
				if (loop == 0) {
					for (int i = 0; i < D; i++) {
						tX[i] = 0;
					}
				}
				else {
					for (int i = 0; i < D; i++) {
						tX[i] = sqrt(th) * rnd.nextGaussian();
					}
				}
				for (int itr = 0; itr < 100 * D; itr++) {
					vector<double> g = grad(tX, 1.0);
					for (int i = 0; i < D; i++) {
						tX[i] -= 0.1 * g[i];
					}

					double dd = dst_from_origin(tX);
					if (dd > th) {
						dd = sqrt(dd);
						for (int j = 0; j < D; j++) {
							tX[j] *= (1 - 1e-2) * (sqrt(th) / dd);
						}
					}
				}
				double score = calc(tX, 1.0);
				if (bestscore > score) {
					bestscore = score;
					best = tX;
				}
			}
			swap(tX, best);

			X = mean;
			for (int i = 0; i < D; i++) {
				for (int j = 0; j < D; j++) {
					X[i] += _ss[i][j] * tX[j];
				}
			}
		}
		else {
			X = Gaussian(mean, sigma, 2.0);
		}
		return X;
	}

	void sample_result(const vector<double> &X, double Y) {

		pq.emplace_back(Y, X);
		int pos = 0;
		for (int i = pq.size() - 1; i > 0; i--) {
			if (pq[i].first < pq[i - 1].first) {
				swap(pq[i], pq[i - 1]);
			}
			else {
				pos = i;
				break;
			}
		}
		curS[NEvals - 1] = pq[0].first;
		while (pq.size() > 200) pq.pop_back();

	}
};

class GP_CMA_ES {
private:

	int NEvals = 0;
	vector<double> X;
	vector<double> mean;
	vector<vector<double> > sigma;
	vector<pair<double, vector<double> > > pq;
	double rate = 1;
	int sz = 5;
	double c1 = 0.1;
	Limited_GP gp;
	int cur = -1;

public:
	
	GP_CMA_ES() {
		X.resize(D, 0);
		init_mean_sigma(mean, sigma);
	}

	GP_CMA_ES(int _sz) {
		sz = _sz;
		init_mean_sigma(mean, sigma);
	}

	pair<double, vector<double> >  get_best_point() {
		return pq[0];
	}

	vector<double> get_next_sample() {

		NEvals++;
		if (cur == -1) {
			cur = 0;
			vector<double> i_mean = Gaussian(mean, sigma, 0.01);
			gp.init(i_mean, sigma);
		}
		vector<double> X = gp.get_next_sample();
		return X;
	}

	void sample_result(const vector<double> &X, double Y) {
		double prog = (double)cnt / N;
		gp.sample_result(X, Y);

		if (gp.get_NEvals() > 10 * (D)) {
			if (gp.converge() || gp.get_NEvals() > 1000 * D) {
				cur = -1;
			}
		}
		if (cur == 0)return;

		pq.push_back(gp.get_best_point());
		int pos = 0;
		for (int i = pq.size() - 1; i > 0; i--) {
			if (pq[i].first < pq[i - 1].first) {
				swap(pq[i], pq[i - 1]);
			}
			else {
				pos = i;
				break;
			}
		}
		if (pos == 0) {
			gp.update_bestS();
		}
		for (int i = pq.size() - 2; i >= 0; i--) {
			if (pq[i + 1].first - pq[i].first < 1e-7) {
				if (dst(pq[i + 1].second, pq[i].second) < 1e-4) {
					pq.erase(pq.begin() + (i + 1));
				}
			}
		}

		while (pq.size() > 200) pq.pop_back();

		if ((int)pq.size() >= 10) {

			vector<vector<double> > v;
			sz = min(10, (int)pq.size());
			for (int i = 0; i < sz; i++) {
				v.push_back(pq[i].second);
			}
			vector<double> W(v.size(), 1);
			double r = pow(1e-9 + 1 - prog, 1.0 / sz);
			for (int i = 0; i < v.size(); i++) {
				W[i] = pow(r, i);
			}
			vector<double> nmean(D);
			vector<vector<double> > nsigma;
			calc_sigma(v, W, nmean, nsigma, false);
			for (int i = 0; i < D; i++) {
				mean[i] = nmean[i];
			}
			for (int i = 0; i < D; i++) {
				for (int j = 0; j < D; j++) {
					sigma[i][j] = (1 - c1) * sigma[i][j] + c1 * nsigma[i][j];
				}
			}
		}
	}

};

int main() {
	input();
	vector<double> X;
	double Y;
	GP_CMA_ES gp_cma_es;
	while (cnt < N) {
		X = gp_cma_es.get_next_sample();
		Y = ask(X);
		gp_cma_es.sample_result(X, Y);
	}
	return 0;
}