#ifndef DNV_INCLUDED_EXTERNIO
#define DNV_INCLUDED_EXTERNIO 1
#include "commons.h"
#include <sys/socket.h>
#include <arpa/inet.h>
#include <unistd.h>
#ifndef NOJSON
#include "json/json.h"
#endif
#ifndef NOZMQ
#include <zmq.hpp>
#include <typeinfo>
#endif
#ifndef PORT
#define PORT 5555
#endif
#ifndef IPADDR
#define IPADDR "192.168.4.25"
#endif

/**This exception is thrown when loading charges externally if the data does not match the format specified.*/
class WrongFormatSpecifiedException : public exception {};
#ifndef NOJSON
/** This exception is thrown if a particular atom is not present in JSON data*/
class AtomNotInJSONException: public exception {}; //When Atom type variable is not initialized in an Atom object
/**
  This is the jsonio namespace. It contains functions for processing JSON data by using the jsoncpp library.<br/>
  In case DeNovo is compied without JSON support, this namespace (along with its functions is not available for use)

  <br/><i>jsoncpp is a standard library to process data in JSON format.</i> See https://en.wikibooks.org/wiki/JsonCpp for a brief introduction to the library.<br/>https://github.com/open-source-parsers/jsoncpp is the GitHub repository.
*/
namespace jsonio
{
  /**@brief Read a JSON format string to produce a usable Json::Value object.*/
  Json::Value& readJSON(const std::string& is)
  {
    Json::CharReaderBuilder builder;
    Json::CharReader * reader = builder.newCharReader();
    static std::string errors;
    Json::Value* ret=new Json::Value();
    if(!reader->parse(is.c_str(),is.c_str()+is.length(),ret,&errors)) cout << errors<<"\n";
    return *ret;
  }
  /**@brief Same as readJSON(const std::string&) but uses a stream for input instead of a string*/
  inline Json::Value& readJSON(istream& is) {return readJSON(stringfx::drain(is));}
  /**@brief Write JSON data to screen*/
  void dumpJSON(Json::Value& v)
  {
    for(const auto& el : v)
      cout << el << "\n";
    //for(int i=0;i<v.size();i++) cout << v[(int)i] << "\n";
  }
}
#endif
/**
  This is the txtio namespace. This namespace has functions to process text in different formats (comma-separated, and JSON string) to load into standard vectors which can be understood by DeNovo (See: Molecule::assignChargesFromData())
*/
namespace txtio
{
  /**@brief Load charges from a comma-separated text input. For the expected input format see algorithm PDF*/
  std::vector<std::pair<std::string,double>> loadChargesFromText(const std::string& ins)
  {
    if(ins.find('}')!=std::string::npos) throw WrongFormatSpecifiedException();
    std::vector<std::pair<std::string,double>> ret;
    std::vector<std::string> lines=stringfx::split(ins,',');
    std::string n; double v;
    for(std::string& s : lines)
    {
      s=stringfx::trim(s);
      stringstream ss(s);
      ss >> n >> v;
      ret.push_back(make_pair(n,v));
    }
    return ret;
  }
  /**@brief Load charges from a JSON String. For the expected input format see algorithm PDF*/
  std::vector<std::pair<std::string,double>> loadChargesFromJSONString(const std::string& jstr)
  {
    std::vector<std::pair<std::string,double>> ret;
    int curInd=jstr.find('{')+1;
    if(curInd==std::string::npos) throw WrongFormatSpecifiedException();
    int endInd=jstr.find('}');
    curInd=jstr.find('\"',curInd);
    int tempI;
    std::string atid; double chg;
    do
    {
      curInd++;
      tempI=curInd;
      curInd=jstr.find('\"',curInd);
      if(jstr[curInd-1]=='\\') atid=jstr.substr(tempI,curInd-tempI-1);
      else atid=jstr.substr(tempI,curInd-tempI);
      curInd=jstr.find(':',curInd)+1;
      tempI=curInd;
      curInd=jstr.find(',',curInd);
      if(curInd==std::string::npos) curInd=jstr.find('\\',tempI);
      chg=std::stod(jstr.substr(tempI,curInd-tempI));
      ret.push_back(make_pair(atid,chg));
      curInd=jstr.find('\"',curInd+1);
      if(curInd>=endInd) break;
    }
    while(curInd!=std::string::npos);
    return ret;
  }
}
/**
  The netio namespace contains functions to perform I/O operations over a network.
  It uses native c++ sockets for communication.<br/>
  See: zmqio (the namespace)
*/
namespace netio
{
  static int sock = 0;
  /**@brief Prepare a socket for communication. There can only be a single communication socket in a given run.*/
  int prepareSocket(const std::string& ip,std::ostream& dout=std::cout)
  {
    dout << "Prepare socket called\n";
    struct sockaddr_in serv_addr;
    if ((sock = socket(AF_INET, SOCK_STREAM, 0)) < 0)
    {
        cout << "ERR: Socket creation error \n";
        return 0;
    }

    serv_addr.sin_family = AF_INET;
    serv_addr.sin_port = htons(PORT);

    // Convert IPv4 and IPv6 addresses from text to binary form
    if(inet_pton(AF_INET, ip.c_str(), &serv_addr.sin_addr)<=0)
    {
        cout << "ERR: Invalid address/ Address not supported \n";
        return 0;
    }

    if (connect(sock, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0)
    {
        cout << "ERR: Connection Failed \n";
        return 0;
    }
    dout << "Connected to socket\n";
    return sock;
  }
  /**@brief Send a message over socket (as a string) an receive the reply*/
  std::string communicateWithSocket(int sock,const std::string& inp,std::ostream& dout=std::cout)
  {
    char buffer[1024] = {0};
    send(sock , inp.c_str() , inp.length() , 0 );
    //cout << "Sent. Awaiting response\n";
    int valread = recv( sock , buffer, 1024,0);
    dout << buffer << "\n";
    //cout << "Response received\n";
    return std::string(buffer);
  }
  /**@brief Get the initialized socket*/
  inline int getSocket() {return sock;}
  /**@brief Initialize a new socket with the given IPADDR and PORT. These can be defined using "#define"s*/
  inline void init(std::ostream& dout=std::cout) {prepareSocket(IPADDR,dout);}
}
#ifndef NOZMQ
/**
  This is the zmqio namespace. It contains methods similar to the netio namespace, but uses the ZMQ library for communication.<br/>
  When compiled without ZMQ support, this namespace (with all iths functions) is disabled. Using ZMQ for communication assumes that the server also uses ZMQ for I/O over the network.

  <br/><i>This namespace uses the ZMQ communication libraries</i>. See: https://zeromq.org/languages/cplusplus
*/
namespace zmqio
{
  class ZMQIOThread
  {
    zmq::context_t context{ 1 };
    zmq::socket_t socket;

  public:
    /**@brief Initialize a new socket with the given IPADDR and PORT. These can be defined using "#define"s*/
    ZMQIOThread(std::ostream& dout=std::cout) {this->prepareSocket(IPADDR,"tcp",PORT,dout);}

    /**@brief Get the initialized ZMQ I/O thread object's socket*/
    inline zmq::socket_t& getZMQSocket() {return socket;}
    /**@brief Get the initialized ZMQ I/O thread object's socket (Alternate function for ZMQIOThread::getZMQSocket())*/
    inline zmq::socket_t& getSocket() {return socket;}
    /**@brief Prepare a socket for communication. There can only be a single communication socket in a given object ZMQIOThread object.*/
    zmq::socket_t& prepareSocket(const std::string& ip,const std::string& protocol="tcp",int port=PORT,std::ostream& dout=std::cout)
    {
      dout << "Prepare socket (ZMQ) called\n";
      // construct a REQ (request) socket and connect to interface
      socket={ context, zmq::socket_type::req };
      std::string socket_name=protocol+"://"+ip+":"+to_string(port);
      dout << socket_name << "\n";
      socket.connect(socket_name);
      dout << "Connection established to: "<<ip<<":"<<port<<"\n";
      return socket;
    }
    /**@brief Send a message over socket (as a string) an receive the reply*/
    std::string communicateWithSocket(zmq::socket_t& socket,const std::string& data,std::ostream& dout=std::cout)
    {
      // send the request message
      //std::cout << "Sending Request for charges" << request_num << "..." << std::endl;
      socket.send(zmq::buffer(data), zmq::send_flags::none);
      // wait for reply from server
      zmq::message_t reply{};
      socket.recv(reply, zmq::recv_flags::none);
      std::string re((char*)reply.data(),reply.size());
      dout << re << "\n";
      return re;
    }
  };
  // initialize the zmq context with a single IO thread
  //static zmq::context_t context{ 1 };
  //static zmq::socket_t socket;
}
#endif
#endif
