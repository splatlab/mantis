#include <cstdlib>
#include <iostream>
#include <memory>
#include <utility>
#include <boost/asio.hpp>

#include "query.hh"

using boost::asio::ip::tcp;

class session
  : public std::enable_shared_from_this<session>
{
public:
  session(tcp::socket socket)
    : socket_(std::move(socket))
  {
  }

  void start()
  {
    do_read();
  }

private:
  void do_read()
  {
    // read query file path and output file path
    auto self(shared_from_this());
    socket_.async_read_some(boost::asio::buffer(data_, max_length),
        [this, self](boost::system::error_code ec, std::size_t length)
        {
          std::cout << "Read: " << data_ << std::endl;

          char * delimeter_ptr = strstr(data_, " "); 
          if (!delimeter_ptr) 
          {
            // TODO incorrect query format error. Call do_read?
            return; 
          }
          std::string query_filepath(data_, delimeter_ptr - data_);
          std::string output_filepath(delimeter_ptr + 1);  // assume exactly one space between filepaths

          query::run_query(query_filepath, output_filepath); // TODO check result

          if (!ec)
          {
            do_write(length, query_filepath);
          }
        });
  }

  void do_write(std::size_t length, std::string outfile)
  {
    // respond with output filepath 
    auto self(shared_from_this());
    std::cout << "Respond with message: " << data_ << std::endl;
    boost::asio::async_write(socket_, boost::asio::buffer(data_, length),
        [this, self](boost::system::error_code ec, std::size_t /*length*/)
        {
          if (!ec)
          {
            do_read();
          }
        });
  }

  tcp::socket socket_;
  enum { max_length = 4096 };
  char data_[max_length];
};

class server
{
public:
  server(boost::asio::io_service& io_service, short port)
    : acceptor_(io_service, tcp::endpoint(tcp::v4(), port)),
      socket_(io_service)
  {
    do_accept();
  }

private:
  void do_accept()
  {
    acceptor_.async_accept(socket_,
        [this](boost::system::error_code ec)
        {
          if (!ec)
          {
            std::make_shared<session>(std::move(socket_))->start();
          }

          do_accept();
        });
  }

  tcp::acceptor acceptor_;
  tcp::socket socket_;
};