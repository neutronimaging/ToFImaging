#ifndef FINDCLOSEST_H
#define FINDCLOSEST_H

#include <QDebug>
#include <iostream>
using namespace std;

namespace ToFImagingAlgorithms {


template<typename T>
T getClosest(T, T, T);

// Returns element closest to target in arr[]
template<typename T>
int findClosest(T *arr, int n, T target)
{

    std::cout << arr[0] << std::endl;
    std::cout << target << std::endl;
    std::cout << n << std::endl;
    // Corner cases
    if (target <= arr[0])
    {
        std::cout << " I am returning zero" << std::endl;
        return 0;
    }
    if (target >= arr[n - 1])
        return n-1;

    // Doing binary search
    int i = 0, j = n, mid = 0;

    std:: cout << "before while " << std::endl;
    while (i < j) {
        mid = (i + j) / 2;
        std::cout<< "mid: " << mid << std::endl;

        if (arr[mid] == target)
        {
            std::cout << " I am returning mid: " << mid << std::endl;
            return mid;
        }

        /* If target is less than array element,
            then search in left */
        if (target < arr[mid]) {

            // If target is greater than previous
            // to mid, return closest of two
            if (mid > 0 && target > arr[mid - 1])
                return getClosest(arr[mid - 1],
                                  arr[mid], target);

            /* Repeat for left half */
            j = mid;

            std::cout << "j: " <<j << std::endl;
        }

        // If target is greater than mid
        else {
            if (mid < n - 1 && target < arr[mid + 1])
                return getClosest(arr[mid],
                                  arr[mid + 1], target);
            // update i
            i = mid + 1;
            std::cout << "i: " << i << std::endl;
        }
    }

    std:: cout << "after while " << std::endl;
    std::cout << mid << std::endl;

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
T getClosest(T val1, T val2, T target)
{
    if (target - val1 >= val2 - target)
        return val2;
    else
        return val1;
}

}
#endif // FINDCLOSEST_H
