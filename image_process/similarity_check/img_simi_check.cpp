/*
author:Franc_zi
time_stamp:2020_04_21
description:This is a beta version without test.You could modify and input the initial data.
refferd_docs:图形相似原理_谭建荣
*/
#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#define TOP 2
#define GEO 4
#define SCALE 5
using namespace std;

typedef struct item
{
	int id;
	int link;
	int scale;
}Item;

typedef struct image
{
	int view, id, aux;
	vector<Item> seq;
}Img;

void classify_img(vector<Img>& seq, unordered_map<int, vector<Img>>& list) {
	for (auto& it : seq) list[it.id].push_back(it);
}

double element_wise_type(Img& left, Img& right) {
	auto seq_l = left.seq;
	auto seq_r = right.seq;
	double total = 0;
	static unordered_set<int> s1 = { 0,1,2,3 };
	for (int i = 0; i < seq_l.size(); i++) {
		auto l_id = seq_l[i].id; auto r_id = seq_r[i].id;
		if (l_id == r_id) total += 1;
		else if (s1.count(l_id) && s1.count(r_id)) total += 0.25;
	}
	return total;
}

double linked_state(Img& left, Img& right) {
	auto seq_l = left.seq;
	auto seq_r = right.seq;
	int size = seq_l.size();
	bool l_judge, r_judge;
	double total = 0;
	for (int i = 0; i < size; i++) {
		l_judge = r_judge = false;
		auto l_id_next = seq_l[(i + 1) % size].id;
		auto r_id_next = seq_r[(i + 1) % size].id;
		if (l_id_next == r_id_next) total += 1;
	}
	return total;
}

double linked_type(Img& left, Img& right) {
	auto seq_l = left.seq;
	auto seq_r = right.seq;
	double total = 0;
	static unordered_set<int> s1 = { 0,6 }, s2 = { 1,2,5 };
	for (int i = 0; i < seq_l.size(); i++) {
		auto l_link = seq_l[i].link; auto r_link = seq_r[i].link;
		if (l_link == r_link) total += 1;
		else if (s1.count(l_link) && s2.count(r_link)) total += 0.75;
		else if ((s1.count(l_link) && s2.count(r_link)) || (s1.count(r_link) && s2.count(l_link))) total += 0.5;
	}
	return total;
}

double symmetry_check(Img& left, Img& right) {
	if (left.aux == right.aux) return 1;
	else return 0;
}

double scale_check(Img& left, Img& right) {
	auto seq_l = left.seq;
	auto seq_r = right.seq;
	double total = 0;
	static unordered_set<int> s1 = { 2,3 };
	for (int i = 0; i < seq_l.size(); i++) {
		auto l_scale = seq_l[i].scale; auto r_scale = seq_r[i].scale;
		if (l_scale == r_scale) total += 1;
		else if (s1.count(l_scale) && s1.count(r_scale)) total += 0.75;
	}
	return total;
}

double cal_similarity(int code, Img& item1, Img& item2) {
	int left, right;
	double up = 0, down = 1;
	if (code == TOP) down--;
	else if (code == GEO) down += 3 * item1.seq.size();
	else down += 4 * item1.seq.size();
	for (int i = 0; i < code; i++) {
		if (i == 0) up += element_wise_type(item1, item2);
		else if (i == 1) up += linked_state(item1, item2);
		else if (i == 2) up += linked_type(item1, item2);
		else if (i == 3) up += symmetry_check(item1, item2);
		else if (i == 4) up += scale_check(item1, item2);
	}
	return up / down;
}

int main()
{
	//scan the images and format to Img for acquiring the seq manually;
	vector<Img> seq;
	unordered_map<int, vector<Img>> list;
	classify_img(seq, list);
	for (auto& it : list) {
		auto seq = it.second;
		auto size = seq.size();
		for (int i = 0; i < size - 1; i++) {
			for (int j = i + 1; j < size; j++) {
				cout << "(" << i << " " << j << " " << "):" << cal_similarity(SCALE, seq[i], seq[j]) << endl;
				cout << endl;
			}
		}
		cout << endl;
	}
	return 0;
}
