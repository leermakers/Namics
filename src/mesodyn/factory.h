
#ifndef FACTORY_H
#define FACTORY_H

#include <memory>
#include <map>

template< class Base, typename Key_type, class... Args >
class Factory_template
{
public:

   typedef std::shared_ptr<Base> (*Factory_function_type)(Args...);

   Factory_template(const Factory_template&) = delete;
   Factory_template &operator=(const Factory_template&) = delete;

   static void Register(const Key_type &key, Factory_function_type fn)
   {
      function_list[key] = fn;
   }

   static std::shared_ptr<Base> Create(const Key_type &key, Args... args)
   {
      auto iter = function_list.find(key);
      if (iter == function_list.end())
         return 0;
      else
         return (iter->second)(args...);
   }

   static Factory_template &Instance() { static Factory_template gf; return gf; }

private:
   Factory_template() = default;

   typedef std::map<Key_type, Factory_function_type> Function_map;
   static Function_map function_list;
};

template< class Base, typename Key_type, class... Args >
typename Factory_template<Base, Key_type, Args...>::Function_map Factory_template<Base, Key_type, Args...>::function_list;

template <class Base, class Derived, typename Key_type, class... Args>
class Register_class
{
public:
   Register_class(const Key_type &key)
   {
      Factory_template<Base, Key_type, Args...>::Register(key, factory_function);
   }

   static std::shared_ptr<Base> factory_function(Args... args)
   {
      return std::make_shared<Derived>(args...);
   }
};

#endif