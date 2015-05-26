//
//  MergeList.h
//  
//
//  Created by Giovanni on 11/06/2013.
//
//

#ifndef _MergeList_h
#define _MergeList_h

#include <list>




template <class T>
class MergeList: public std::list<T> {
	
	//Sub-class of std::list::iterator that implements
	//also next() and previous() functions
	template <class Q>
	class NextIterator : public std::list<Q>::iterator {
		
		typedef typename std::list<T>::iterator iterator_base;
		
	public:
		
		//Default constructor
		NextIterator<Q> () : iterator_base() {};
		
		//Copy constructor
		NextIterator<Q> (const iterator_base& r) :
		iterator_base(r) {};
		
		//Assignment operator
		NextIterator<Q>& operator = (const iterator_base& rhs)
		{
			iterator_base::operator =(rhs);
			return *this;
		}
		
		//Next
		NextIterator<Q> next() {
			NextIterator<Q> temp = *this;
			temp++;
			return temp;
		}
		
		//Previous
		NextIterator<Q> previous() {
			NextIterator<Q> temp = *this;
			temp--;
			return temp;
		}
		
	};
	
	
	typedef typename std::list<T> list_base;
public:

	typedef NextIterator<T> iterator;
    
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
	
	iterator begin() {
		iterator it = list_base::begin();
		return it;
	}
	
	iterator end() {
		iterator it = list_base::end();
		return it;
	}
	
};

#endif
