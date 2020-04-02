#include <iostream>
#include <cfloat>
#include <cmath>
using namespace std;
double r1,r2,a,b,c,d,e,f,tmp;

typedef struct coordinates {
	double x;
	double y;
	double z;
}VT;

inline double power(double base) {return base*base;}

inline double linear_quadratic_fc(double a,double b,double c,double x) {return a*x*x+b*x+c;}

inline bool judge(double x,double lbd) {return x>0&&x<1&&lbd>=0;}

double kt_valid(double x) {
	static int cnt = 0;
	double lbd,res = DBL_MAX;
	if(cnt==0) {
		lbd = -f*x - d;
		if(judge(x,lbd)) res = linear_quadratic_fc(c,e,a,x);
	}else if(cnt==1) {
		lbd = e - f*x;
		if(judge(x,lbd)) res = linear_quadratic_fc(b,-d,a,x);
	}else if(cnt==2) {
		lbd = -(2*c+e-f*x);
		if(judge(x,lbd)) res = linear_quadratic_fc(b,-d-f,a+c+e,x);
	}else if(cnt==3) {
		lbd = d-2*b+f*x;
		if(judge(x,lbd)) res = linear_quadratic_fc(c,e-f,a+b-d,x);
	}
	cnt++;
	return res;
}

double cal(int id) {
	switch(id) {
		case 0: return a;
		case 1: return a+c+e;
		case 2: return a+b-d;
		case 3: return a+b+c-d+e-f;
		case 4: return kt_valid(-e/(2*c));
		case 5: return kt_valid(d/(2*b));
		case 6: return kt_valid((d+f)/(2*b));
		case 7: return kt_valid((f-e)/(2*c));
		case 8: 
			double x2 = (f*d-2*b*e)/(4*b*c-f*f);
			double x1 = (d+f*x2)/(2*b);
			if(x1>0&&x1<1&&x2>0&&x2<1) return a + b*x1*x1 + c*x2*x2 - d*x1 + e*x2 - f*x1*x2;
			else return DBL_MAX;
	}
}

int main() {
	double min_val=DBL_MAX;
	VT v[4];
	cout << "please enter the two radius respectively:" << endl;
	cin >> r1 >> r2;
	cout << "please enter the four coordinates(i.e.for a,b,c,d)" << endl;
	for(int i=0;i<4;i++) cin >> v[i].x >> v[i].y >> v[i].z;
	a = power(v[0].x-v[2].x) + power(v[0].y-v[2].y) + power(v[0].z-v[2].z);
	b = power(v[1].x-v[0].x) + power(v[1].y-v[0].y) + power(v[1].z-v[0].z);
	c = power(v[3].x-v[2].x) + power(v[3].y-v[2].y) + power(v[3].z-v[2].z);
	d = 2 * (v[2].x-v[0].x) * (v[1].x-v[0].x) +
		2 * (v[2].y-v[0].y) * (v[1].y-v[0].y) +
		2 * (v[2].z-v[0].z) * (v[1].z-v[0].z);
	e = 2 * (v[2].x-v[0].x) * (v[3].x-v[2].x) +
		2 * (v[2].y-v[0].y) * (v[3].y-v[2].y) +
		2 * (v[2].z-v[0].z) * (v[3].z-v[2].z);
	f = 2 * (v[1].x-v[0].x) * (v[3].x-v[2].x) +
		2 * (v[1].y-v[0].y) * (v[3].y-v[2].y) +
		2 * (v[1].z-v[0].z) * (v[3].z-v[2].z);
	for(int i=0;i<9;i++) min_val = min(min_val,((tmp=cal(i))>=0?tmp:DBL_MAX));
	min_val = sqrt(min_val);
	if(min_val>=r1+r2) 
		cout << "There will be safe between these two robot arms!" << endl;
	else 
		cout << "There will be a crash between these two robot arms!" << endl;
	cout << "The min dist val: " << min_val << endl; 
	return 0;
}
