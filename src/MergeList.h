//
//  MergeList.h
//  
//
//  Created by Giovanni on 11/06/2013.
//
//

#ifndef _MergeList_h
#define _MergeList_h

#include<list>


template <class T>
class MergeList: public std::list<T> {
    
public:

    typedef typename std::list<T>::iterator iterator;
    
    //This function overloads the operator + between two lists,
    //so that it's possible to concatenate them (as it happens
    //in Python)
    MergeList operator+ (MergeList<T>& list2) {
        MergeList<T> result;
        iterator it;
        
        //Add elements from this
        for (it = this->begin(); it != this->end(); it++) {
            result.push_back(*it);
        }
        
        //Add elements from the other list
        for (it = list2.begin(); it != list2.end(); it++) {
            result.push_back(*it);
        }
        
        return result;
    }
};

#endif
