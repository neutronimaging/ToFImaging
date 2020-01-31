#ifndef FINDCLOSEST_H
#define FINDCLOSEST_H

#include <iostream>


namespace ToFImagingAlgorithms {


template<typename T>
int getClosest(T, T, T);

// Returns element closest to target in arr
template<typename T>
unsigned int findClosest(T *arr, int n, T target)
{

    // Corner cases
    if (target <= arr[0])
    {
        return 0;
    }
    if (target >= arr[n - 1])
        return n-1;

    // Doing binary search
    int i = 0, j = n, mid = 0;

    while (i < j) {
        mid = (i + j) / 2;

        if (arr[mid] == target)
            {
            return mid;
        }

        /* If target is less than array element,
            then search in left */
        if (target < arr[mid]) {

            // If target is greater than previous
            // to mid, return closest of two
            if (mid > 0 && target > arr[mid - 1])
            {

               if (getClosest(arr[mid - 1], arr[mid], target)==0)
                    {
                    return mid-1;
                }
                else {
                    return mid;
                }

            }

            /* Repeat for left half */
            j = mid;

        }

        // If target is greater than mid
        else {
            if (mid < n - 1 && target < arr[mid + 1])
            {
                if (getClosest(arr[mid], arr[mid + 1], target)==0)
                    return mid;
                else {
                    return mid+1;
                }
            }
            // update i
            i = mid + 1;
        }
    }


    // Only single element left after search
    // returns the array position
    return mid;
}

// Method to compare which one is the more close.
// We find the closest by taking the difference
// between the target and both values. It assumes
// that val2 is greater than val1 and target lies
// between these two.
template<typename T>
unsigned int getClosest(T val1, T val2, T target)
{
    if (target - val1 >= val2 - target)
        return 1;
    else
        return 0;
}

}
#endif // FINDCLOSEST_H
