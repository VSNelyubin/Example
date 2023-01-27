
#include <iostream>
#include <vector>
#include <fstream>

typedef double db;
typedef std::vector<db> da;
typedef std::vector<da> dm;

db spv(db a, int x) {
	db rez = 1;
	for (int i = 0; i < x; i++) {
		rez *= a;
	}
	return rez;
}

db fz(db x, db y, db z) {
	return 1 / x * z + 2 / x * y;
}

db fzz(da& x, da& y, da& z,int i) {
	return 1 / x[i] * z[i] + 2 / x[i] * y[i];
}

db f(db x) {
	return (cos(2) - sin(2)) * cos(2 * sqrt(x)) + (cos(2) + sin(2)) * sin(2 * sqrt(x));
}

db ff(db x) {
	return sin(x);
}

db f5(db x) {
	return 0;
}

db inter(db x0,db x1,db s5, db (*func)(db)) {
	db rez = 0;
	for (db x5 = x0 + s5; x5 < x1; x5 += s5 * 2) {
		rez += s5 / 3 * (func(x5 - s5) + func(x5 + s5) + 4 * func(x5));
		//std::cout << "added " << x5 << ", " << x5 + s5 << " and " << x5 + s5 / 2 << "\n";

	}
	return rez;
}

double U(double x, double y, double z) {
	return -4 * x * z - (4 * x * x + 2) * y;
}

double V(double z) {
	return z;
}

double f2(double x) {
	return (1 + x) * exp(-x * x);
}

double p(double x) {
	return -4 * x;
}

double q(double x) {
	return -(4 * x * x + 2);
}

		//y'' +1/x y' +2/x y=0
		//y(1)=1
		//y'(1)=1
		//x[1,2] h-0.1

		//z' -1/x * z-2/x y=0
		//z'=1/x*z + 2/x y
		//z(1)=1
		//y(1)=1
//=============================================
		//y''+4xy'+(4x^2+2)y=0
		//y'=z
		
		//z'+4xz+(4x^2+2)y=0
		//z'=-4xz-(
		//y'(0)=1
		//4y(2)-z(2)=23e^(-4)
		// 
		//

da progon(dm cofas) {
	int lpl = cofas[0].size();
	double tp = 0;
	da p;
	da q;
	da x;
	x.push_back(0);
	x.push_back(0);
	tp = -cofas[2][0] / cofas[1][0];
	p.push_back(tp);
	tp = cofas[3][0] / cofas[1][0];
	q.push_back(tp);
	for (int i = 1; i < lpl; i++) {
		tp = -cofas[2][i] / (cofas[1][i] + cofas[0][i] * p[i - 1]);
		p.push_back(tp);
		tp = (cofas[3][i] - cofas[0][i] * q[i - 1]) / (cofas[1][i] + cofas[0][i] * p[i - 1]);
		q.push_back(tp);
		x.push_back(0);
	}


	for (int j = 0; j < lpl; j++) {
		int i = lpl - 1 - j;
		x[i] = p[i] * x[i + 1] + q[i];
	}
	x.pop_back();
	return x;
}




int main() {
	std::cout.precision(5);
	int launcher = 3;
	if (launcher & 1) {
		std::cout << "part 1\n";
		//y'=z
		std::ofstream myfile;
		myfile.open("example.txt");
		da y1;
		da z1;
		z1.push_back(1);
		y1.push_back(1);
		db h1 = 0.1;
		//int i = 0;
		db dl;
		da x;
		std::cout << "Euler:\n";
		int si = 1.0 / h1;
		for (int i = 0; i < si; i++) {
			x.push_back(1 + h1 * i / si);
			dl = z1[i];
			y1.push_back(y1[i] + dl * h1);
			dl = -fz(x[i], y1[i], z1[i]);
			z1.push_back(z1[i] + dl * h1);
			std::cout << y1[i] << "\t err = " << f(1 + i * h1) - y1[i] << "\n";
			//myfile << "ar["<<i-1<<"]=" << y1[i] << "\n";
		}


		da z2, y2;
		z2.push_back(1);
		y2.push_back(1);
		db q1, q2, q3, q0;
		std::cout << "\nRunge-Khutt:\n";
		db k1, k2, k3, k0;
		for (int i = 0; i < si; i++) {
			k0 = z2[i];
			q0 = fz(x[i], y2[i], z2[i]);
			k1 = z2[i] + q0 * h1 / 2;
			q1 = fz(x[i] + h1 / 2, y2[i] + k0 * h1 / 2, z2[i] + q0 * h1 / 2);
			k2 = z2[i] + q1 * h1 / 2;
			q2 = fz(x[i] + h1 / 2, y2[i] + k1 * h1 / 2, z2[i] + q1 * h1 / 2);
			k3 = z2[i] + q2 * h1;
			q3 = fz(x[i], y2[i] + k2 * h1, z2[i] + q2 * h1);


			dl = h1 / 6 * (k0 + 2 * k1 + 2 * k2 + k3);
			y2.push_back(y2[i] + dl);
			dl = h1 / 6 * (k0 + 2 * k1 + 2 * k2 + k3);
			//dl = -x;
			z2.push_back(z2[i] + dl);
			std::cout << y2[i] << "\t err = " << f(1 + i * h1) - y2[i] << "\n";
			//myfile << "ar[" << i - 1 << "]=" << y1[i] << "\n";
		}


		da z3, y3;
		std::cout << "\nAdams:\n";

		for (int ii = 0; ii < 4; ii++) {
			z3.push_back(z2[ii]);
			y3.push_back(y2[ii]);
		}
		for (int i = 0; i < si; i++) {
			if (z3.size() > i) {
				//std::cout << z3[i] << "\n";// << f(x) << "\n";
				//myfile << "arr[" << i << "]=" << z3[i] << "\n";
				//i++;
				continue;
			}
			i--;
			dl = h1 * h1 / 360 * (323 * fzz(x, y3, z3, i) - 264 * fzz(x, y3, z3, i - 1) + 159 * fzz(x, y3, z3, i - 2) - 38 * fzz(x, y3, z3, i - 3));

			y3.push_back(y3[i] + dl);
			dl = h1 / 24 * (55 * fzz(x, y3, z3, i) - 59 * fzz(x, y3, z3, i - 1) + 37 * fzz(x, y3, z3, i - 2) - 9 * fzz(x, y3, z3, i - 3));
			//std::cout << dl << "\n";
			//dl = -x;
			z3.push_back(z3[i] + dl);
			i++;
			//i++;
			//std::cout << z3[i] << "\n";// << f(x) << "\n";
			///myfile << "arr[" << i  << "]=" << z3[i] << "\n";
			//myfile << "ar[" << i - 1 << "]=" << y1[i] << "\n";
		}
		for (int i = 3; i < z3.size() - 1; i++) {
			z3[i + 1] = z3[i] - h1 / 720 * (-251 * fzz(x, y3, z3, i + 1) - 646 * fzz(x, y3, z3, i) + 264 * fzz(x, y3, z3, i - 1) + 159 * fzz(x, y3, z3, i - 2) - 38 * fzz(x, y3, z3, i - 3));
			y3[i + 1] = y3[i] - h1 * h1 / 1440 * (-135 * fzz(x, y3, z3, i + 1) - 752 * fzz(x, y3, z3, i) + 246 * fzz(x, y3, z3, i - 1) - 96 * fzz(x, y3, z3, i - 2) + 17 * fzz(x, y3, z3, i - 3));
		}

		for (int i = 0; i < y3.size(); i++) {
			std::cout << y3[i] << "\t err = " << f(1 + i * h1) - y3[i] << "\n";
			//myfile << "arr[" << i  << "]=" << y3[i] << "\n";
		}

		myfile.close();
	}
	if (launcher & 2) {
		std::cout << "part 2\n";
		if (true) {
			//y'=z
			//z'+4xz+(4x^2+2)y=0
			// 
			//y'(0)=1
			//4y(2)-z(2)=23e^(-4)
			std::cout << "metod strelbi\n";
			double nu = 0;
			//da = 1;
			da z2;
			da y2;
			double dl = 0;
			da x2;
			double dlte = 1;
			bool loop = true;
			double h1 = 0;
			double h2 = 2;
			double h = 0.1;
			int si = 1.0 / h;
			int count = 0;
			double paster = 0;
			double newer = 0;
			while (loop) {
				//nu += 0.2;
				z2.clear();
				y2.clear();
				x2.clear();
				z2.push_back(nu);
				y2.push_back(1);

				for (int i = 0; i < si + 1; i++) {
					x2.push_back(h1 + (h2-h1) * i / si);
					double x = x2[i];
					double y = y2[i];
					double z = z2[i];

					double k0 = V(z);

					double q0 = U(x, y, z);

					double k1 = V(z + q0 * h / 2);

					double q1 = U(x + h / 2, y + k0 * h / 2, z + q0 * h / 2);

					double k2 = V(z + q1 * h / 2);

					double q2 = U(x + h / 2, y + k1 * h / 2, z + q1 * h / 2);

					double k3 = V(z + q2 * h);

					double q3 = U(x + h, y2[i] + k2 * h, z + q2 * h);

					dl = h / 6 * (q0 + 2 * q1 + 2 * q2 + q3);
					z2.push_back(z + dl);
					dl = h / 6 * (k0 + 2 * k1 + 2 * k2 + k3);
					y2.push_back(y + dl);

					//print(str(h1 + i * hi)[3:] + ":" + str(y2[i])[3:] + " err = " + str(f2(h1 + i * hi) - y2[i])[3:])
				}


				//y = y2
				//x = np.linspace(0, 2, len(y))

				//count += 0.1

				//color = (math.sin(count) * *2, 0.5 * math.sin(count) * *2, 1 - math.sin(count) * *2)
				//plt.plot(x, y, c = color)
				double ch = 4 * y2[y2.size() - 1] - z2[z2.size() - 1];
				//std::cout << dlte << "\t" << ch - 23 * exp(-4)<<"\n";
				if (abs(ch - 23 * exp(-4)) < 0.01) {
					break;
				}

				newer = abs(ch - 23 * exp(-4));
				if (newer > paster) {
					dlte /= -2;
				}
				nu += dlte;
				paster = newer;
				//plt.plot(x, y, c = color)
				//for i in range(0, len(x)) :
				//	y[i] = f2(x[i])
				//	plt.plot(x, y, c = 'red')
				//	plt.show()
			}
			x2.push_back(h2);
			for (int i = 0; i < y2.size(); i++) {

				std::cout << y2[i] << "\t err = " << f2(x2[i]) - y2[i] << "\n";
			}

		}
		if (true) {
			std::cout << "\n";
			std::cout << "Finite difference method\n";
			da x2;
			da y2;
			db	h1 = 0;
			db h2 = 2;
			db h = 0.2;
			int N = int(1.0 / h);
			da a;
			da b;
			da c;
			da d;
			db ya = 1;
			db yb = 0;
			db x;
			for (int i = 0; i < N + 1; i++) {
				x2.push_back(h1 + (h2-h1) * i / N);
				x = x2[i];
				if (i == 0) {
					a.push_back(0);
					c.push_back(1 + p(x) * h / 2);
					d.push_back(h* h* f(x) - (1 - p(x) * h / 2) * ya);
				}
				else if (i == N) {
					yb = 23 * exp(-4);
					c.push_back(0);
					a.push_back(1 - p(x) * h / 2);
					d.push_back(h* h* f(x) - (1 + p(x) * h / 2) * yb);
				}
				else {
					a.push_back(1 - p(x) * h / 2);
					c.push_back(1 + p(x) * h / 2);
					d.push_back(h* h* f(x));
				}
				b.push_back(-2 + h * h * q(x));
			}
			dm args = { a, b, c, d };

			da y = progon(args);

			yb = (23 * exp(-4) - 1 / h * y[N - 1]) / (4 - 1 / h);
			d[N] = h * h * f(x) - (1 + p(x) * h / 2) * yb;
			y = progon(args);
			for (int i = 0; i < y.size(); i++) {
				std::cout << y[i] << "\t err = " << f2(x2[i]) - y[i] << "\n";
			}
			std::cout << "\n";
		}
	}
	return 0;
}
