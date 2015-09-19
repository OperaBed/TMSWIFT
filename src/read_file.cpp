void read_ofile(std::string file_name, std::vector<double> *vec)
{
        std::string line;
        std::fstream infile;
        infile.open(file_name.c_str());
        while (std::getline(infile, line))
        {
                std::istringstream iss(line);
                double a;
                if ((iss >> a))
                {
                        vec->push_back(a);
                }
        }
        infile.close();
}

