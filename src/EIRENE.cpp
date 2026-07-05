#include "EIRENE.h"
#include <cstdlib>
#define e 2.71828

EIRENE::EIRENE(int Fit, int Num, const string &database)
{
    Fit_ = Fit;
    Num_ = Num;
    database_ = database;
    string line;
    int b;
    string database_file;
    if (database_ == "amjuel")
        database_file = "amjuel.txt"; /*amjuel.tex*/ // path and name of v_ionisation file
    else if (database_ == "methane")
        database_file = "METHANE";
    else if (database == "hydhel")
        database_file = "hydhel.txt"; /*hydhel.tex*/ // path and name of v_ionisation file
    else if (database == "ammonx_elsa")
        database_file = "AMMONX_Arrh-elast.tex";
    string path_EIRENE;
    std::ifstream fp;
    const char *env_database_path = std::getenv("NT_DATABASE_PATH");
    vector<string> database_paths;
    if (env_database_path)
        database_paths.push_back(env_database_path);
    database_paths.push_back("Inputfile/database/");
    database_paths.push_back("nt/Inputfile/database/");
    database_paths.push_back("/data/leuven/379/vsc37950/nt/Inputfile/database/");
    for (const string &database_path : database_paths)
    {
        string candidate = database_path;
        if (!candidate.empty() && candidate.back() != '/')
            candidate += "/";
        candidate += database_file;
        fp.clear();
        fp.open(candidate, std::ios::in);
        if (fp.is_open())
        {
            path_EIRENE = candidate;
            break;
        }
    }

    if (!fp.is_open())
    {
        cout << "This file READING for " + path_EIRENE + " have some problem!!!\n";
    }
    if (Fit_ == 2 || Fit_ == 8)
    {
        string ex;
        int n = 2;
        for (int i = 0; i < Num; i++)
            getline(fp, line);
        // getline(fp, ex);
        /*fp.seekg(n, ios::cur);
        fp >> data[0];
        cout << data[0] << endl;*/
        /*
        for (int i = 0; i < 3; i++) {
            data[i + 3 * j] = ex[5];
        }*/
        for (int i = 0; i < 9; i++)
        {
            fp >> ex;
            if (ex == "\end{verbatim}\end{small}")
                break;
            fp >> data_[i];
            // std::cout << data_[i] << endl;
        }
        /*for (int i = 0; i < 9; i++)
            cout << i << "\t" << data_[i] << endl;*/
    }
    if (Fit_ == 3 || Fit_ == 4 || Fit_ == 10)
    {
        seekline(fp, Num);
        for (int k = 0; k < 3; k++)
        {
            for (int i = 0; i < 9; i++)
            {
                fp >> b;
                for (int y = 0; y < 3; y++)
                {
                    fp >> data_[i * 9 + (k * 3 + y)];
                }
            }
            Num += 12;
            seekline(fp, Num);
        }
    }
    fp.close(); // �ر���
    // out();
}

double EIRENE::cal(double n, double T, double *H2dataforH8)
{
    double v_EIRENE = 0;
    if (Fit_ == 2)
    {
        for (int i = 0; i < 9; i++)
        {
            v_EIRENE += data_[i] * pow(log(T), i);
        }
        v_EIRENE = exp(v_EIRENE);
    }
    if (Fit_ == 3)
    {
        if (n < 0.1)
        {
            n = 0.1;
        }
        if (T < 0.1)
        {
            T = 0.1;
        }
        for (int i = 0; i < 9; i++)
        {
            for (int j = 0; j < 9; j++)
            {
                v_EIRENE = v_EIRENE + data_[9 * i + j] * pow(log(n), j) * pow(log(T), i);
            }
        }
        // std::cout << "EIRENE check: " << this << " " << n << ", " << T << " " << v_EIRENE << ", " << exp(v_EIRENE) << endl;
        v_EIRENE = exp(v_EIRENE);
    }
    if (Fit_ == 4)
    {
        if (T < 0.1)
        {
            T = 0.1;
        }
        n = n / 1.e6;
        if (n <= 1.e8)
        {
            for (int i = 0; i < 9; i++)
            {
                v_EIRENE = v_EIRENE + data_[i * 9] * pow(log(T), i);
            }
        }
        else
        {
            n /= 1.e8;
            for (int i = 0; i < 9; i++)
            {
                for (int j = 0; j < 9; j++)
                {
                    v_EIRENE = v_EIRENE + data_[9 * i + j] * pow(log(n), j) * pow(log(T), i);
                }
            }
        }
        v_EIRENE = exp(v_EIRENE);
    }
    if (Fit_ == 8)
    {
        /*n = n / 1.e6;
        double n3 = 0., V = 3. / 2;
        V = cal(n, T, H2dataforH8);
        for (int i = 0; i < 8; i++)
        {
            n3 += ;
        }
        v_EIRENE = exp(v_EIRENE);*/
        v_EIRENE = 0.896 * T * (*H2dataforH8);
        return v_EIRENE * 1.6e-19 / 1e6;
    }
    if (Fit_ == 10)
    {
        n = n / 1.e6;
        if (n <= 1.e8)
        {
            for (int i = 0; i < 9; i++)
            {
                v_EIRENE = v_EIRENE + data_[i * 9] * pow(log(T), i);
            }
        }
        else
        {
            n /= 1.e8;
            for (int i = 0; i < 9; i++)
            {
                for (int j = 0; j < 9; j++)
                {
                    v_EIRENE = v_EIRENE + data_[9 * i + j] * pow(log(n), j) * pow(log(T), i);
                }
            }
        }
        v_EIRENE = exp(v_EIRENE);
    }
    if (v_EIRENE < 0)
        return 0;
    else
    {
        if (Fit_ == 0)
            return v_EIRENE / 1.e4;
        else
            return v_EIRENE / 1.e6;
    }
}

// C-R file reading.cpp
ifstream &seekline(ifstream &in, int line) // ���򿪵��ļ�in����λ��line�С�
{
    char buf[1024];
    in.seekg(0, ios::beg); // ��λ���ļ���ʼ��
    for (int i = 0; i < line; i++)
    {
        in.getline(buf, sizeof(buf)); // ��ȡ�С�
    }
    return in;
    in.close();
}

double EIRENE::data(int i)
{
    return data_[i];
}

int EIRENE::fit() const
{
    return Fit_;
}

void EIRENE::out()
{
    cout << Fit_ << '\t' << Num_ << '\t';
    if (Fit_ == 2 || Fit_ == 8)
    {
        for (int i = 0; i < 9; i++)
            cout << data_[i] << '\t';
    }
    if (Fit_ == 3 || Fit_ == 4 || Fit_ == 10)
    {
        for (int i = 0; i < 81; i++)
            cout << data_[i] << '\t';
    }
    cout << endl;
}
