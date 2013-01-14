#include "tools.h"
#include <vector>
using namespace std;

// TODO - change to work with new data alignment!
void quickSort(hitData *items, int left, int right)
{
  int i, j;
  double x;
  hitData tmp;

  i = left;
  j = right;
  x = items[(left+right)/2].x;

  do {
    while((items[i].x < x) && (i < right))
       i++;
    while((x < items[j].x) && (j > left))
       j--;

    if(i <= j) {
    	//trackHitId zamienić
      tmp = items[i];
      items[i] = items[j];
      items[i].inputHit = i;//
      items[j] = tmp;
      items[j].inputHit = j;//
      i++; j--;
    }
  } while(i <= j);

  if(i < right)
      quickSort(items, i, right);
  if(left < j)
      quickSort(items, left, j);
}

void sortHits(hitData *hits, int sensNum, sensorInfo sensors[]){
	int start = 0;
	for(int i=0; i<sensNum; i++){
		quickSort(hits,start,start+sensors[i].hitsNum-1);
		start+=sensors[i].hitsNum;
	}
}


