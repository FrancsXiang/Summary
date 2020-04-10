/*
	author:Franc_zi
	time_stamp:2020_04_10
	refer_docs:https://blog.csdn.net/wangjoe11/article/details/78019578
*/
#include <iostream>
#include <cstdio>
using namespace std;

template<typename T>
T* getContainer(int len) {new T[len];}

void THA(double* m,double* h,double* y,int n,bool natural=true) {
	double* r = getContainer<double>(n+2);
	double* a = getContainer<double>(n+2);
	double* b = getContainer<double>(n+2);
	double* c = getContainer<double>(n+2);
	double* g = getContainer<double>(n+2);
	double* d = getContainer<double>(n+2);
	for(int i=2;i<=n;i++) {
		r[i] = 6 * ((y[i]-y[i-1])/h[i-1]-(y[i-1]-y[i-2])/h[i-2]);
		a[i] = h[i-2];
		b[i] = 2 * (h[i-2]+h[i-1]);
		c[i] = h[i-1];
	}
	if(natural) {
		b[1] = b[n+1] = 1;
		c[1] = a[n+1] = 0;
		r[1] = r[n+1] = 0;
	}else {
		cout << "please enter the start and end first order derivative:" << endl;
		double A,B;
		cin >> A >> B;
		b[1] = 2 * h[0];c[1] = h[0];
		a[n+1] = h[n-1];b[n+1] = 2 * h[n-1];
		r[1] = 6 * ((y[1]-y[0])/h[0] - A);
		r[n+1] = 6 * (B - (y[n]-y[n-1])/h[n-1]); 
	}
	g[1] = c[1] / b[1];
	d[1] = r[1] / b[1];
	for(int i=2;i<=n;i++) {
		if(i<n) g[i] = c[i] / (b[i] - g[i-1] * a[i]);
		d[i] = (r[i] - d[i-1] * a[i]) / (b[i] - g[i-1] * a[i]);
	}
	m[n] = g[n+1];
	for(int i=n-1;i>=0;i--) m[i] = d[i+1] - g[i+1] * m[i+1];
}

int main()
{
	int n;
	double a,b,c,d;
	cout << "please enter the total number of shape points:" << endl;
	cin >> n;
	cout << "please enter the coordinates respectively:(i.e.(x,y))" << endl;
	double* x = getContainer<double>(n);
	double* y = getContainer<double>(n);
	for(int i=0;i<n;i++) cin >> x[i] >> y[i];
	n--;
	double* h = getContainer<double>(n);
	double* m = getContainer<double>(n+1);
	for(int i=0;i<n;i++) h[i] = x[i+1] - x[i];
	THA(m,h,y,n);
	for(int i=0;i<n;i++) {
		a = y[i];
		b = (y[i+1]-y[i])/h[i] - h[i]*m[i]/2 - h[i]*(m[i+1]-m[i])/6;
		c = m[i] / 2;
		d = (m[i+1]-m[i]) / (6*h[i]);
		printf("The %d interval params:(%.2f,%.2f,%.2f,%.2f)\n",i+1,a,b,c,d);
	}
	return 0;
}
