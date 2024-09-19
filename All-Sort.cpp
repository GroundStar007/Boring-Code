#include <iostream>
#include <cstdio>
#include <algorithm>
#include <vector>
#include <stack>
#include <cstring>

using namespace std;

// 插入排序  O(n^2)  稳定 
inline void InsertSort(int *a, int n) {
	for (int i = 1; i < n; i++) {
		int end = i;
		int x = a[end + 1];
		while (end >= 1) {
			if (a[end] > x) {
				a[end + 1] = a[end];
				end --;
			}
			else
				break;
		}
		a[end + 1] = x;
	}
}

// 希尔排序（优化插入）  O(n^2)  稳定 
inline void ShellSort(int *a, int n) {
	int gap = n / 3 + 1;
	for (int i = 1; i <= n - gap; i++) {
		int end = i;
		int x = a[end + gap];
		while (end >= 1) {
			if (a[end] > x) {
				a[end + gap] = x;
				end -= gap;
			}
			else break;
		}
		a[end + gap] = x;
	} 
} 

// 选择排序  O(n^2)  不稳定 
inline void SelectSort(int *a, int n) {  // 与数据无关，任何数据都是 O(N^2) 
	int begin = 1;
	int end = n;
	while (begin < end) {
		int maxi = a[begin];
		int mini = a[begin];
		for (int i = begin; i <= end; i++) {
			if (a[i] > a[maxi])
				maxi = i;
			if (a[i] < a[mini])
				mini = i;
		}
		swap(a[begin], a[mini]);
		if (maxi == begin)
			maxi = mini;
		swap(a[maxi], a[end]);
		begin ++;
		end --;
	}
}

// 堆排序  O(nlogn) 不稳定  
//（也许涉及交换的算法都不稳定） 
inline void AdjustDown(int* a, int k, int root) { 
	int parent = root;
	int child = parent * 2;
	while (child <= k) {
		if (child + 1 <= k && a[child + 1] > a[child])
			child += 1;
		if (a[child] > a[parent]) {
			swap(a[child], a[parent]);
			parent = child;
			child = parent * 2;
		}
		else break;
	}
} 
 
inline void HeapSort(int *a, int n) {
	for (int i = n / 2 - 1; i >= 1; i--) 
		AdjustDown(a, n, i);
	for (int i = n; i > 1; i--) {
		swap(a[1], a[i]);
		AdjustDown(a, i - 1, 1);  // 注意为 i-1，为最后一个没有被排序的位置 
	}
}

// 冒泡排序  O(n^2) 稳定
// 稳定是因为只有大于才会交换，等于不会交换 
inline void BubbleSort(int *a, int n) {
	bool flag = false;
	for (int i = 1; i < n; i++) {
		for (int j = 1; j < n - i + 1; j++) {  // 此处小于是因为需要使用 a[j+1] 
			if (a[j] > a[j + 1]) {
				// swap 内部装入元素值 
				swap(a[j], a[j + 1]);
				flag = true;  // 记录发生交换，数组还未有序 
			}
		}
		// 如果没有发生过交换，数组已经有序，不必再执行 
		if (flag == false)
			break; 
	}
} 

// hoare 方法进行排序 
inline int PartSort1(int* a, int left, int right) {
	int key = right;
	while (left < right) {
		while (left < right && a[left] <= a[key]) left ++;
		while (left < right && a[right] >= a[key]) right --;
		swap(a[left], a[right]);	
	}
	swap(a[left], a[key]);
	return left;
}

// 挖坑法 
inline int PartSort2(int* a, int left, int right) {
	int key = a[left];
	int hole = left;
	while (left < right) {   // 先更新右侧，再更新左侧 
		while (left < right && a[right] >= key) right --;
		a[hole] = a[right];
		hole = right;

		while (left < right && a[left] <= key) left ++;
		a[hole] = a[left];
		hole = left;
	}
	a[hole] = key;
	return hole;
}

// 双指针法
inline int PartSort3(int* a, int left, int right) {
	int key = left;
	int cur = left + 1;
	int prev = left;
	while (cur <= right) {
		if (a[cur] <= a[key] && ++prev < cur)
			swap(a[prev], a[cur]);
		++ cur;
	}
	swap(a[prev], a[key]);
	return prev;
}

// 三数取中法取基准值（为了防止取到最小值而产生不必要的循环） 
inline int MidiumIndex(int* a, int left, int right) {
	vector<int> v;
	v.push_back(a[left]);
	v.push_back(a[right]);
	v.push_back(a[(left + right) >> 1]);
	sort(v.begin(), v.end());
	return v[1];
}

// 快速排序（递归方法） 
// O(nlogn) 对于基本有序的序列最劣，退化到 O(n^2)
// 不稳定  
// 在较小区间内可尝试改为插入排序以降低时间复杂度 
inline void QuickSort(int* a, int left, int right) {
	if (left > right) return;
	int key = PartSort3(a, left, right);
	QuickSort(a, left, key - 1);
	QuickSort(a, key + 1, right);
}

// 非递归法实现快速排序 
// !! 要注意栈后进先出的特性，先压入右端点 
inline void QuickSortNonR(int* a, int left, int right) {
	stack<int> s;
	s.push(right);
	s.push(left);
	while (!s.empty()) {
		int L = s.top(); s.pop();
		int R = s.top(); s.pop();
		int mid = PartSort1(a, L, R);
		if (R > mid + 1) {  // 判定是否符合要求，隐性判定 
			s.push(R);
			s.push(mid + 1);
		}
		if (L < mid - 1) {
			s.push(mid - 1);
			s.push(L);
		}
	}
}

// 归并排序 递归实现 
// O(nlogn)  稳定 
inline void MergeSort(int* a, int left, int right) {
	if (left >= right) 
		return;
		
	int mid = (left + right) >> 1;
	MergeSort(a, left, mid);		//左区间
	MergeSort(a, mid + 1, right);	//右区间
	
	int b[100];  // 需要再使用一个数组 
	int cnt = left;    // 计数器 
	int begin1 = left, end1 = mid;
	int begin2 = mid + 1, end2 = right;
	while (begin1 <= end1 && begin2 <= end2) {
		if (a[begin1] < a[begin2])  // 注意使左端点自加 
			b[cnt++] = a[begin1++];
		else
			b[cnt++] = a[begin2++]; 
	}
	while (begin1 <= end1) b[cnt++] = a[begin1++];
	while (begin2 <= end2) b[cnt++] = a[begin2++];
	for (int i = left; i <= right; i++)
		a[i] = b[i];
}

// 非递归实现归并排序 
inline void MergeSortNonR(int* a, int n) {
	int b[100];
	int gap = 1;
	while (gap < n) {
//		memset(b, 0, sizeof b);
		//归并每两组归并一次
		int index = 0; //记录b数组中的元素下标
		for (int i = 1; i <= n; i += 2 * gap) {	//两组中的元素个数为 2 * gap
			
			int begin1 = i, end1 = i + gap - 1;	//边界
			int begin2 = i + gap, end2 = i + 2 * gap - 1;
			
            //当原数组中元素个数不是 2^n 时，最后两组组会出现元素不匹配的情况
            //情况1：end1 >= n 或 begin2 >= n,即最后两组中只有一组有元素，则不需归并
			if (end1 > n || begin2 > n) break;
            //情况2：end2 >= n，即最后两组中，第二组元素个数小于第一组，则需要调整第二组边界
			if (end2 > n) end2 = n;
			while (begin1 <= end1 && begin2 <= end2) {
				if (a[begin1] < a[begin2])
					b[++index] = a[begin1++];
				else
					b[++index] = a[begin2++];
			}
			while (begin1 <= end1) b[++index] = a[begin1++];
			while (begin2 <= end2) b[++index] = a[begin2++];
		}
        //一趟排完后，将归并后的有序序列拷贝到原数组中
	    for (int j = 1; j <= index; j++) a[j] = b[j];
		gap *= 2;
		// 真的要注意输入的坐标从 1 开始还是从 0 开始啊，太麻烦了 
	}
}
 
const int N = 1e5 + 10;
int n, a[N];

int main() {
	cin >> n;
	for (int i = 1; i <= n; i++)
		cin >> a[i];
		
//	InsertSort(a, n);
//	ShellSort(a, n);
//	SelectSort(a, n);
//	HeapSort(a, n);
//	BubbleSort(a, n);
//	QuickSort(a, 1, n);
//	QuickSortNonR(a, 1, n);
//	MergeSort(a, 1, n);
	MergeSortNonR(a, n); 
	
	for (int i = 1; i <= n; i++)
		cout << a[i] << " ";
	return 0;
}