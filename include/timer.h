/*====================================================
|
|
|¼ÇÊ±Æ÷
|
|
=====================================================*/
#ifndef _TIMER_H
#define _TIMER_H


#include <time.h>
#include <assert.h>
#include <iostream>

class Timer {

public:
    Timer() 
    { 
        state_ = uninitialized;
    }

    void start()
    { 
        state_ = running;
        t1_ = systemTime();
    }

    void stop()
    {
        t2_ = systemTime();
        assert(state_ == running);
        state_ = stopped;
    }   
    long double elapsedSeconds()
    {
        assert(state_ == stopped);
        return t2_ - t1_;
    }
	
	void diplayelapsedtime()
	{
		std::cout<<elapsedSeconds()<<"s\n";
	}
private:
    Timer(Timer&) { }

    long double systemTime()
    {
        return clock() / (long double) CLOCKS_PER_SEC;
    }

    enum { uninitialized, running, stopped } state_;

    long double t1_, t2_;
};



#endif //_TIMER_H

