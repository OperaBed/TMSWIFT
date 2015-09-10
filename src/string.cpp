std::string toStr(unsigned long long num) {
   std::string charset = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz+-.,?<>;:_=!@#$^&*()~`";
   int base = charset.length();
   std::string str = num ? "" : "0";

   while (num) {
      str = charset.substr(num % base, 1) + str;
      num /= base;
   }

   return str;
}

