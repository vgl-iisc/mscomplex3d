#ifndef __UTL_H_INCLUDED
#define __UTL_H_INCLUDED

/*===========================================================================*/
/* Misc utility functions
/*---------------------------------------------------------------------------*/

#include <vector>
#include <string>
#include <sstream>

namespace utl {

/*---------------------------------------------------------------------------*/

///  \brief A generic to string impl that uses the << operator of type T
template <typename T>
inline std::string to_string (const T& v)
{std::stringstream ss; ss << v ; return ss.str();}

/*---------------------------------------------------------------------------*/

template <typename T1,typename T2>
inline std::string to_string (const std::pair<T1,T2> &p)
{std::stringstream ss; ss <<p.first <<","<< p.second ; return ss.str();}

/*---------------------------------------------------------------------------*/

template <class iter_t>
void argsort(iter_t b, iter_t e, std::vector<size_t>& idxs);

/*---------------------------------------------------------------------------*/

template <class iter_t>
void rearrange(iter_t b, iter_t e,const std::vector<size_t>& idxs);

/*---------------------------------------------------------------------------*/

/// \brief strip front and back whitespaces in given string
void trim(std::string &s);

/*---------------------------------------------------------------------------*/

/// \brief Converts a given string to the template type
template <class T> inline T from_string(std::string s);

/*---------------------------------------------------------------------------*/

///  \brief Converts a given string to the template types
template <typename T1,typename T2>
inline void from_string(std::string s,T1& t1,T2& t2)
{std::stringstream ss(s);ss >> t1 >> t2;}

/*---------------------------------------------------------------------------*/

///  \brief Converts a given string to the template types
template <typename T1,typename T2,typename T3>
inline void from_string(std::string s,T1& t1,T2& t2,T3& t3)
{std::stringstream ss(s);ss >> t1 >> t2 >> t3;}

/*---------------------------------------------------------------------------*/

///  \brief Converts a given string to the template types
template <typename T1,typename T2,typename T3,typename T4>
inline void from_string(std::string s,T1& t1,T2& t2,T3& t3,T4& t4)
{std::stringstream ss(s);ss >> t1 >> t2 >> t3 >> t4;}

/*---------------------------------------------------------------------------*/

}// namespace utl
/*===========================================================================*/






/*===========================================================================*/
/* Misc utility classes
/*---------------------------------------------------------------------------*/

//#include <boost/iterator/iterator_facade.hpp>
//#include <boost/shared_ptr.hpp>
//#include <boost/date_time/posix_time/posix_time.hpp>
//#include <boost/thread/mutex.hpp>
//#include <boost/thread/thread.hpp>

#include <iostream>
#include <string>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <mutex>

#ifdef __linux__
#include<memory>
#endif

namespace utl {

/*---------------------------------------------------------------------------*/


/**
  \brief A class to iterate over non-comment lines of a file

  \note All parts of a line that succeed the comment char ('#' default) are
  assumed to be a comment and are stripped

  \note Sample usage
  for(file_line_iterator lgen("file.txt"),lend; lgen != lend; )
    cout << *lgen++;

**/

class file_line_iterator
    /* : public boost::iterator_facade<
    file_line_iterator
    , std::string const
    , boost::forward_traversal_tag
    >*/
{
public:
    

  file_line_iterator(const char * f,char c_char='#');
	file_line_iterator();

private:
  //friend class boost::iterator_core_access;
  void increment();
  bool equal(file_line_iterator const& other) const;
  const std::string &dereference() const;

private:
  //boost::shared_ptr<std::ifstream> is;
  std::shared_ptr<std::ifstream> is;
  std::string                      value;
  char                             c_char;
};

	/*---------------------------------------------------------------------------*/

/** \brief Simple stopwatch timer (to measure wall clock time NOT cpu time) **/

	class timer
{
 public:
  timer()
  {restart();}

  inline void   restart()
  {
	//  _start_time = boost::posix_time::microsec_clock::local_time();
      _start_time = std::chrono::steady_clock::now();

  }

  inline double elapsed() const
  {
    //boost::posix_time::time_duration td =
      //  boost::posix_time::microsec_clock::local_time() - _start_time;

    auto now = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_time = now - _start_time;
    return elapsed_time.count();
    //return double(td.total_milliseconds())/1000;
  }

 private:
  //boost::posix_time::ptime _start_time;
  std::chrono::steady_clock::time_point _start_time;

}; // timer


/*---------------------------------------------------------------------------*/

/** \brief A simple multi-threaded logger                                   **/

class logger {
public:
    enum severity_level { trace, debug, info, warning, error, fatal };
	inline bool isOpen(severity_level severity) {
        return severity >= UTL_LOG_LEVEL;
    }

    inline void push_ts(const std::string& log) {
        auto tstr = get_current_time();
        std::lock_guard<std::mutex> lock(s_mutex);
        std::clog << "[" << tstr << "] " << log << std::endl;
    }

    inline void push(const std::string& log) {
        std::lock_guard<std::mutex> lock(s_mutex);
        std::clog << log;
    }

    static inline logger& get() {
        return s_logger;
    }

private:
    static std::string get_current_time() {
        auto now = std::chrono::system_clock::now();
        auto in_time_t = std::chrono::system_clock::to_time_t(now);
        auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
            now.time_since_epoch()) % 1000;

        std::ostringstream oss;
        oss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %H:%M:%S");
        oss << '.' << std::setfill('0') << std::setw(3) << milliseconds.count();
        return oss.str();
    }

    static std::mutex s_mutex;
    static logger s_logger;
};




namespace detail
{
/** \brief a small extension of std::pair where only first is initialized   **/
template<class _T1, class _T2> struct pair:public std::pair<_T1,_T2>
{pair(const _T1 &t1){std::pair<_T1,_T2>::first = t1;}};
}//namespace detail

}// namespace utl

/*===========================================================================*/



/*===========================================================================*/
/* Binary read and write functions
/*---------------------------------------------------------------------------*/

namespace utl {

/*---------------------------------------------------------------------------*/

template<typename T>
inline void bin_write(std::ostream &os, const T &v)
{os.write((const char*)(const void*)&v,sizeof(T));}

/*---------------------------------------------------------------------------*/

template<typename T>
inline void bin_read(std::istream &is, const T &v)
{is.read((char*)(void*)&v,sizeof(T));}

/*---------------------------------------------------------------------------*/

template<typename T>
inline void bin_read_raw(std::istream &is, T * v, size_t n)
{is.read((char*)(void*)v,n*sizeof(T));}

/*---------------------------------------------------------------------------*/

template<typename T>
inline void bin_write_raw(std::ostream &os, T * v, size_t n)
{os.write((const char*)(const void*)v,n*sizeof(T));}

/*---------------------------------------------------------------------------*/

template<typename T>
inline void bin_read_vec(std::istream &is, std::vector<T> &v)
{
  size_t n;
  bin_read(is,n);
  v.resize(n);
  is.read((char*)(void*)v.data(),v.size()*sizeof(T));
}

/*---------------------------------------------------------------------------*/

template<typename T>
inline void bin_write_vec(std::ostream &os, const std::vector<T> &v)
{
  bin_write<size_t>(os,v.size());
  os.write((const char*)(const void*)v.data(),v.size()*sizeof(T));
}

/*---------------------------------------------------------------------------*/

}
/*===========================================================================*/


/*===========================================================================*/
/* Macro definitions
/*---------------------------------------------------------------------------*/

#define is_in_range(i,b,e) (((b) <= (i)) && ((i) < (e)))
#define two_power(i)   (std::pow(2,(i)))
#define divide_rz(a,b) ((a)/(b))
#define divide_ri(a,b) (((a)+(b)-1)/(b))

/*---------------------------------------------------------------------------*/

#define LOG(lev)  \
  if(utl::logger::get().isOpen(utl::logger::lev))\
  for(utl::detail::pair<bool,std::stringstream> __utl_lm_v__(true); \
      __utl_lm_v__.first ;__utl_lm_v__.first=false,\
  utl::logger::get().push(__utl_lm_v__.second.str())) __utl_lm_v__.second



/*--------------------*/

// Define a macro for different compilers
#ifdef _MSC_VER
#define MY_PRETTY_FUNCTION __FUNCSIG__  // MSVC uses __FUNCSIG__
#elif defined(__GNUC__) || defined(__clang__)
#define MY_PRETTY_FUNCTION __PRETTY_FUNCTION__  // GCC and Clang use __PRETTY_FUNCTION__
#else
#define MY_PRETTY_FUNCTION __func__  // Fallback: __func__ is part of C++ standard
#endif

/*---------------------------------------------------------------------------*/

#define ELOG  \
  if(utl::logger::get().isOpen(utl::logger::error))\
  for(utl::detail::pair<bool,std::stringstream> __utl_lm_v__(true); \
      __utl_lm_v__.first ;__utl_lm_v__.first=false,\
  utl::logger::get().push_ts(__utl_lm_v__.second.str())) \
  __utl_lm_v__.second \
  <<utl::detail::__classFunction__(MY_PRETTY_FUNCTION) << " "

/*---------------------------------------------------------------------------*/

#define FLOG  \
  if(utl::logger::get().isOpen(utl::logger::fatal))\
  for(utl::detail::pair<bool,std::stringstream> __utl_lm_v__(true); \
      __utl_lm_v__.first ;__utl_lm_v__.first=false,\
  utl::logger::get().push_ts(__utl_lm_v__.second.str())) \
  __utl_lm_v__.second \
  <<utl::detail::__classFunction__(MY_PRETTY_FUNCTION) << " "



/*---------------------------------------------------------------------------*/

#define ILOG  \
  if(utl::logger::get().isOpen(utl::logger::info))\
  for(utl::detail::pair<bool,std::stringstream> __utl_lm_v__(true); \
      __utl_lm_v__.first ;__utl_lm_v__.first=false,\
  utl::logger::get().push_ts(__utl_lm_v__.second.str())) \
  __utl_lm_v__.second \
  <<utl::detail::__classFunction__(MY_PRETTY_FUNCTION) << " "

/*---------------------------------------------------------------------------*/

#define DLOG  \
  if(utl::logger::get().isOpen(utl::logger::debug))\
  for(utl::detail::pair<bool,std::stringstream> __utl_lm_v__(true); \
      __utl_lm_v__.first ;__utl_lm_v__.first=false,\
  utl::logger::get().push_ts(__utl_lm_v__.second.str())) \
  __utl_lm_v__.second \
  <<utl::detail::__classFunction__(MY_PRETTY_FUNCTION) << " "

/*---------------------------------------------------------------------------*/

#define TLOG  \
  if(utl::logger::get().isOpen(utl::logger::trace))\
  for(utl::detail::pair<bool,std::stringstream> __utl_lm_v__(true); \
      __utl_lm_v__.first ;__utl_lm_v__.first=false,\
  utl::logger::get().push_ts(__utl_lm_v__.second.str())) \
  __utl_lm_v__.second << utl::detail::__trace_indenter_t__::get_indent()\
  <<utl::detail::__classFunction__(MY_PRETTY_FUNCTION) << " "


#define TLOG_INDENT utl::detail::__trace_indenter_t__ __trace_indenter__; \

/*---------------------------------------------------------------------------*/

#define LOG_TS(lev)  \
  if(utl::logger::get().isOpen(utl::logger::lev))\
  for(utl::detail::pair<bool,std::stringstream> __utl_lm_v__(true); \
      __utl_lm_v__.first ;__utl_lm_v__.first=false,\
  utl::logger::get().push_ts(__utl_lm_v__.second.str())) __utl_lm_v__.second

/*---------------------------------------------------------------------------*/

#define SVAR(VAR) #VAR << " = "<< (VAR) << " "

/*---------------------------------------------------------------------------*/

#ifndef ASSERT
#ifndef NDEBUG
#define ASSERT(cond) if (!(cond))\
{ std::stringstream ss; \
  ss<<"Failed to assert condition " << #cond <<"\n" \
    <<"at ("<<__FILE__<<","<<__func__<<","<<__LINE__<<") \n "\
  ;throw std::runtime_error(ss.str());}
#else  //ifndef NDEBUG
#define ASSERT(cond)
#endif // ifndef NDEBUG
#endif // ifndef ASSERT

/*---------------------------------------------------------------------------*/

#ifndef ASSERTV
#ifndef NDEBUG
#define ASSERTV(cond,var1) if (!(cond))\
{ std::stringstream ss; \
  ss<<"Failed to assert condition " << #cond <<"\n" \
    <<"at ("<<__FILE__<<","<<__func__<<","<<__LINE__<<") \n "\
    << #var1 <<" = "<< (var1) << "\n"\
  ;throw std::runtime_error(ss.str());}
#else  //ifndef NDEBUG
#define ASSERTV(cond,var1)
#endif // ifndef NDEBUG
#endif // ifndef ASSERTV

/*---------------------------------------------------------------------------*/

#ifndef ASSERTV2
#ifndef NDEBUG
#define ASSERTV2(cond,var1,var2) if (!(cond))\
{ std::stringstream ss; \
  ss<<"Failed to assert condition " << #cond <<"\n" \
    <<"at ("<<__FILE__<<","<<__func__<<","<<__LINE__<<") \n "\
    << #var1 <<" = "<< (var1) << "\n"\
    << #var2 <<" = "<< (var2) << "\n"\
  ;throw std::runtime_error(ss.str());}
#else  //ifndef NDEBUG
#define ASSERTV2(cond,var1,var2)
#endif // ifndef NDEBUG
#endif // ifndef ASSERTV2

/*---------------------------------------------------------------------------*/

#ifndef ASSERTV3
#ifndef NDEBUG
#define ASSERTV3(cond,var1,var2,var3) if (!(cond))\
{ std::stringstream ss; \
  ss<<"Failed to assert condition " << #cond <<"\n" \
    <<"at ("<<__FILE__<<","<<__func__<<","<<__LINE__<<") \n "\
    << #var1 <<" = "<< (var1) << "\n"\
    << #var2 <<" = "<< (var2) << "\n"\
    << #var3 <<" = "<< (var3) << "\n"\
  ;throw std::runtime_error(ss.str());}
#else  //ifndef NDEBUG
#define ASSERTV3(cond,var1,var2)
#endif // ifndef NDEBUG
#endif // ifndef ASSERTV2

/*---------------------------------------------------------------------------*/

#ifndef ASSERTS
#ifndef NDEBUG
#define ASSERTS(cond)\
  if(!(cond)) \
  for(std::stringstream ss ; true ; throw std::runtime_error(ss.str())) \
  ss<<"Failed to ensure condition " << #cond <<"\n" \
    <<"at ("<<__FILE__<<","<<__func__<<","<<__LINE__<<") \n "
#else  //ifndef NDEBUG
#define ASSERTS(cond)\
if(false) \
std::stringstream()
#endif // ifndef NDEBUG
#endif // ifndef ASSERT


/*---------------------------------------------------------------------------*/

#define ENSURE(cond,mesg) if (!(cond))\
{ std::stringstream ss; \
  ss<<"Failed to ensure condition " << #cond <<"\n" \
    <<"at ("<<__FILE__<<","<<__func__<<","<<__LINE__<<") \n "\
    <<"Message : " << #mesg << "\n"\
  ;throw std::runtime_error(ss.str());}

/*---------------------------------------------------------------------------*/

#define ENSUREV(cond,mesg,var1) if (!(cond))\
{ std::stringstream ss; \
  ss<<"Failed to ensure condition " << #cond <<"\n" \
    <<"at ("<<__FILE__<<","<<__func__<<","<<__LINE__<<") \n "\
    <<"Message : " << #mesg << "\n"\
    << #var1 <<" = "<< (var1) << "\n"\
  ;throw std::runtime_error(ss.str());}

/*---------------------------------------------------------------------------*/

#define ENSUREV2(cond,mesg,var1,var2) if (!(cond))\
{ std::stringstream ss; \
  ss<<"Failed to ensure condition " << #cond <<"\n" \
    <<"at ("<<__FILE__<<","<<__func__<<","<<__LINE__<<") \n "\
    <<"Message : " << #mesg << "\n"\
    << #var1 <<" = "<< (var1) << "\n"\
    << #var2 <<" = "<< (var2) << "\n";\
  ;throw std::runtime_error(ss.str());}

/*---------------------------------------------------------------------------*/

#define ENSUREV3(cond,mesg,var1,var2,var3) if (!(cond))\
{ std::stringstream ss; \
  ss<<"Failed to ensure condition " << #cond <<"\n" \
    <<"at ("<<__FILE__<<","<<__func__<<","<<__LINE__<<") \n "\
    <<"Message : " << #mesg << "\n"\
    << #var1 <<" = "<< (var1) << "\n"\
    << #var2 <<" = "<< (var2) << "\n"\
    << #var3 <<" = "<< (var3) << "\n";\
  ;throw std::runtime_error(ss.str());}

/*---------------------------------------------------------------------------*/

#define ENSURES(cond) \
  if(!(cond)) \
  for(std::stringstream ss ; true ; throw std::runtime_error(ss.str())) \
  ss<<"Failed to ensure condition " << #cond <<"\n" \
    <<"at ("<<__FILE__<<","<<__LINE__<<") \n "

/*---------------------------------------------------------------------------*/

#define BTRACE_ERROR(func) \
try { func ; } catch (const std::exception &e) {\
  if(utl::logger::get().isOpen(utl::logger::error)){\
    std::stringstream ss; ss \
    <<"From: ("<<__FILE__<<","<<__func__<<","<<__LINE__<<") \n ";\
    utl::logger::get().push(ss.str());}\
  throw;}

/*---------------------------------------------------------------------------*/

#define CHECK_ERRORV1(func,var1) \
try { func ; } catch (const std::exception &e) {\
  if(utl::logger::get().isOpen(utl::logger::error)){\
    std::stringstream ss; ss \
    <<"From: ("<<__FILE__<<","<<__func__<<","<<__LINE__<<") \n "\
    << #var1 <<" = "<< (var1) << "\n";\
    utl::logger::get().push(ss.str());}\
  throw;}

/*---------------------------------------------------------------------------*/

#define CHECK_ERRORV2(func,var1,var2) \
try { func ; } catch (const std::exception &e) {\
  if(utl::logger::get().isOpen(utl::logger::error)){\
    std::stringstream ss; ss \
    <<"From: ("<<__FILE__<<","<<__func__<<","<<__LINE__<<") \n "\
    << #var1 <<" = "<< (var1) << "\n"\
    << #var2 <<" = "<< (var2) << "\n";\
    utl::logger::get().push(ss.str());}\
  throw;}

/*---------------------------------------------------------------------------*/

#define CHECK_ERRORV3(func,var1,var2,var3) \
try { func ; } catch (const std::exception &e) {\
  if(utl::logger::get().isOpen(utl::logger::error)){\
    std::stringstream ss; ss \
    <<"From: ("<<__FILE__<<","<<__func__<<","<<__LINE__<<") \n "\
    << #var1 <<" = "<< (var1) << "\n"\
    << #var2 <<" = "<< (var2) << "\n"\
    << #var3 <<" = "<< (var3) << "\n";\
    utl::logger::get().push(ss.str());}\
  throw;}

/*---------------------------------------------------------------------------*/

#define CHECK_ERRORS(func,strm) \
try { func ; } catch (const std::exception &e) {\
  if(utl::logger::get().isOpen(utl::logger::error)){\
    std::stringstream ss; ss \
    <<"From: ("<<__FILE__<<","<<__func__<<","<<__LINE__<<") \n "\
    strm <<std::endl;\
    utl::logger::get().push(ss.str());}\
  throw;}

/*---------------------------------------------------------------------------*/

namespace utl{ namespace detail{

std::string __classFunction__(const std::string& prettyFunction);

class __trace_indenter_t__
{
  static int   s_indent;
public:
  __trace_indenter_t__()   {++s_indent;}
  ~__trace_indenter_t__()  {--s_indent;}
  static std::string get_indent();
};

}}


/*===========================================================================*/


#endif
