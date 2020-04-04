/*
	author:Franc_zi
	timestamp:2020_04_04
	referred_paper:广义hough变换_黎自强(链接: https://pan.baidu.com/s/1aAfC9KLVNK4mknqWws6lng 提取码: 39n4)
	There are some bugs need to be fixed,may be you can help me!!!
*/
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <cmath>
#include <ctime>
#define MAXN 50
#define NOISE 50
using namespace std;
typedef pair<int,int> PII;
// pramameters
int cnt;
double tf,nc,td,tk,tl,te,lbd,mc,mt;
vector<PII> v;
double random_uniform() {
	static long int seed = rand() % 1048576;
	seed = 2045 * seed + 1;
	seed %= 1048576;
	return (double)seed / 1048576.0; 
}

template<typename T>
T* createMatrix(int len) {
	T* res = (T*)malloc(len*len*sizeof(T));
	return res;
}

template<typename T>
T power(T a) {
	return a*a;
}

template<typename T>
void print_matrix(T* img) {
	for(int i=0;i<MAXN;i++) {
		for(int j=0;j<MAXN;j++) {
			if(j) putchar(' ');
			cout << img[i*MAXN+j];
		}
		cout << endl << endl;
	}
}

void make_circle(int x,int y,int radius,int* img) {
	if(x >= MAXN) x = MAXN - 1;
	else if(x < 0) x = 0;
	if(y >= MAXN) y = MAXN - 1;
	else if(y < 0) y = 0;
	if(radius < 0) radius = 1;
	else if(radius > MAXN / 2) radius = (MAXN - 1) / 2;
	if(x < radius || x + radius > MAXN - 1 || y < radius || y + radius > MAXN - 1) {
		cout << "Error:can't draw partly of the circle!" << endl;
		return;
	} 
	double r;
	for(int i=0;i<MAXN;i++) {
		for(int j=0;j<MAXN;j++) {
			r = sqrt(power(i-y)+power(j-x));
			if(r>=radius&&r<=sqrt(2)*radius) img[i*MAXN+j] = 1,v.push_back({i,j});
		}
	} 
}

int* pre_process() {
	int x,y,radius,n_cnt,index;
	n_cnt = 0;
	tf = 500;tk = 1.4;te = lbd = 0.5;tl = 1;td = MAXN - 1;
	int* img = createMatrix<int>(MAXN);
	srand(time(0));
	for(int i=0;i<MAXN;i++) {
		for(int j=0;j<MAXN;j++) {
			if(img[i*MAXN+j] != 1) img[i*MAXN+j] = 0;  
			if(n_cnt < NOISE) {
				index = (int)(random_uniform() * MAXN * MAXN);
				img[index] = 1;//add noise points
				v.push_back({i,j});
				n_cnt++;
			}
		}
	}
	cout << "please enter the three params of circle(i.e.(x,y,radius))" << endl;
	cin >> x >> y >> radius;
	make_circle(x,y,radius,img);
	return img;
}

bool is_noise(PII p,int* img) {
	#ifdef debug
	cout << 1 << endl;
	#endif
	bool flag = true;
	int acc1,acc2;
	int off_set[4][8] = {{2,3,2,3,1,4,1,4},{1,2,3,4,0,1,0,1},
			    {2,3,1,2,1,4,0,3},{1,2,1,2,0,3,0,3}};
	int dx[5] = {-2,-1,0,1,2};
	int dy[5] = {-2,-1,0,1,2};
	for(int i=1;i<4;i++) {
		for(int j=1;j<4;j++) {
			int y = p.first + dy[i];
			int x = p.second + dx[j];
			if(x<0||x>=MAXN||y<0||y>=MAXN) continue;
			if(img[y*MAXN+x] == 1 && (y != p.first || x != p.second)) {
				flag = false;
				break;
			}
		}
	}
	if(flag) return flag;
	for(int k=0;k<4;k++) {
		acc1 = acc2 = 0;
		for(int i=off_set[k][0];i<=off_set[k][1];i++) {
			for(int j=off_set[k][2];j<=off_set[k][3];j++) {
				int y = p.first + dy[i];
				int x = p.second + dx[j];
				if(x<0||x>=MAXN||y<0||y>=MAXN) continue;
				acc1 += img[y*MAXN+x];
			}
		}
		for(int i=off_set[k][4];i<=off_set[k][5];i++) {
			for(int j=off_set[k][6];j<=off_set[k][7];j++) {
				int y = p.first + dy[i];
				int x = p.second + dx[j];
				if(x<0||x>=MAXN||y<0||y>=MAXN) continue;
				acc2 += img[y*MAXN+x];
			}
		}
		if(acc1>1&&acc1<=4&&acc1==acc2) return true;
	}
	return false;
}

inline double coef(int cnt,int i,int j) {
	if(cnt == 0) return j;
	else if(cnt == 1) return j*j;
	else if(cnt == 2) return i;
	else if(cnt == 3) return i*j;
}

bool slope_test(PII p,int* img,int slope) {
	#ifdef debug
	cout << 2 << endl;
	#endif
	double param[4],a;
	memset(param,0,sizeof(param));
	int dx[3] = {-1,0,1};
	int dy[3] = {-1,0,1};
	for(int z=0;z<4;z++) {
		for(int i=0;i<3;i++) {
			for(int j=0;j<3;j++) {
				int y = p.first + dy[i];
				int x = p.second + dx[j];
				if(x<0||x>=MAXN||y<0||y>=MAXN) continue;
				param[z] += coef(z,y,x) * img[y*MAXN+x];
			}
		}
	}
	a = (9*param[3]-param[0]*param[2]) / (9*param[1]-param[0]*param[0]);
	if(fabs(slope-a)>tk) return true;
	else return false;
}
bool sides_split_test(PII p,PII p_m,int k,int* img) {
	#ifdef debug
	cout << 3 << endl;
	#endif
	int l0,l1,l2;
	double direct;
	l0 = l1 = l2 = 0;
	int dx[3] = {-1,0,1};
	int dy[3] = {-1,0,1};
	for(int i=0;i<3;i++) {
		for(int j=0;j<3;j++) {
			int y = p.first + dy[i];
			int x = p.second + dx[j];
			if(x<0||x>=MAXN||y<0||y>=MAXN) continue;
			direct = (y-p.second-k*(x-p.first)) * (p_m.first-p.second-k*(p_m.second-p.first));
			if(direct>0) l1++;
			else if(direct<0) l2++;
			else l0++;
		}
	}
	if(l0+l1-l2<tl) return true;
	else return false;
}

bool edge_test(PII p1,PII p2,PII p3,double& a,double& b,double& r) {
	#ifdef debug
	cout << 4 << endl;
	#endif
	double x1,y1,x2,y2,x3,y3,dist;
	y1 = p1.first;x1 = p2.second;
	y2 = p2.first;x2 = p2.second;
	y3 = p2.first;x3 = p2.second;
	a = ((x2*x2+y2*y2-x1*x1-y1*y1)*(y3-y1)-(x3*x3+y3*y3-x1*x1-y1*y1)*(y2-y1))/
		(2*((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)));
	b = ((x2-x1)*(x3*x3+y3*y3-x1*x1-y1*y1)-(x3-x1)*(x2*x2+y2*y2-x1*x1-y1*y1))/
		(2*((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)));
	r = (sqrt(power(x1-a)+power(y1-b))+sqrt(power(x2-a)+power(y2-b))+sqrt(power(x3-a)+power(y3-b))) / 3;
	double low,high;
	low = min(2*r+2,1.1*r);
	high = max(2*r+2,1.1*r);
	for(int i=0;i<v.size();i++) {
		auto p = v[i];
		dist = power(p.first-b) + power(p.second-a);
		if(dist<low||dist>high||fabs(sqrt(dist)-r)>=te) continue;
		mc++;
	}
	mt = 2 * M_PI * r * lbd;
	if(mc > mt) return false;
	else return true;
}

void single_circle_detect(int* img) {
	PII pi;
	double k1,k2,a,b,r;
	for(int i=0;i<v.size();i++) {
		for(int j=i+1;j<v.size();j++) {
			auto p1 = v[i];
			auto p2 = v[j];
			if(power(v[i].first-v[j].first)+power(v[i].second-v[j].second)>td
			|| v[i].second == v[j].second) continue;
			PII p3{(p1.first+p2.first)/2,(p1.second+p2.second)/2};
			k1 = (v[i].first-v[j].first) / (v[i].second - v[j].second);
			k2 = - 1.0 / k1;
			for(int z=0;z<v.size();z++) {
				cnt++;
				cout << a << " " << b << " " << r << " " << cnt << endl;
				if(v[z] != p3) {
					pi.second = v[z].second;
					pi.first = k2 * (p3.second-v[z].second) + p3.first;
					if(pi.first<0||pi.first>=MAXN||is_noise(pi,img)
					||slope_test(pi,img,k1)||sides_split_test(pi,p3,k1,img)||edge_test(p1,p2,p3,a,b,r)) continue;
					else {
						printf("The equation of circle:(x%c%d)^2+(y%c%d)^2=%d^2",(a>=0?'+':'-'),a,(b>=0?'+':'-'),b,r);
						exit(EXIT_SUCCESS);
					}	
				}
			}
		}
	}
	cout << "There is no circle" << endl;
}

int main()
{
	//create a blank image(1 represent pixel, 0 present nothing)
	int* img = pre_process();
	//print_matrix<int>(img);
	single_circle_detect(img);
	return 0;
}
