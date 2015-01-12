#include "utl.h"

#include <sstream>
#include <fstream>

#include <vector>

#include <tr1/functional>


/*===========================================================================*/

using namespace std;

namespace utl {
/*---------------------------------------------------------------------------*/

void trim(std::string &s)
{
  s.erase(s.begin(), std::find_if
          (s.begin(), s.end(),
           std::not1(std::ptr_fun<int, int>(std::isspace))));
  s.erase(std::find_if
          (s.rbegin(), s.rend(),
           std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
}

/*---------------------------------------------------------------------------*/

file_line_iterator::file_line_iterator(const char * f,char c_char):
  is(new std::ifstream(f)),c_char(c_char)
{
  ENSUREV(is->is_open(),"cannot read the file",f);
  increment();
}

/*---------------------------------------------------------------------------*/

void file_line_iterator::increment()
{
  value.clear();

  while(is && value.size() == 0)
  {
    if(is->eof())
    {
      is.reset();
      value.clear();
      break;
    }

    std::getline(*is,value);

    int p = value.find(c_char);

    if ( p < value.size() )
      value = value.substr(0,p);

    trim(value);
  }
}

/*---------------------------------------------------------------------------*/

bool file_line_iterator::equal(file_line_iterator const& other) const
{
  if(!is && !other.is)
    return true;

  if(!is || !other.is)
    return false;

  ENSURE(is == other.is, "cannot compare distinct istream iters");

  return is->tellg() == other.is->tellg();
}

/*---------------------------------------------------------------------------*/

const std::string & file_line_iterator::dereference() const
{
  ENSURE(is,"dereferencing past end of line stream");
  return value;
}

/*---------------------------------------------------------------------------*/

namespace detail{
std::string __classFunction__(const std::string& prettyFunction)
{
  std::string str(prettyFunction);

  str = str.substr(str.find(" ")+1);
  str = str.substr(0,str.find("("));

  size_t first_colon  =  str.rfind("::");

  if(first_colon != std::string::npos)
  {
    size_t second_colon = str.substr(0,first_colon).rfind("::");

    if(second_colon != std::string::npos)
      str   = str.substr(second_colon+2);
  }
  return str;
}

/*---------------------------------------------------------------------------*/

std::string __trace_indenter_t__::get_indent()
{
  std::string s;

  for(int i = 0 ; i < s_indent; ++i)
    s += "  ";

  return s;
}


}

/*---------------------------------------------------------------------------*/

boost::mutex logger::s_mutex;
logger       logger::s_logger;

int    detail::__trace_indenter_t__::s_indent = 0;

/*---------------------------------------------------------------------------*/

} // namespace utl

/*===========================================================================*/

