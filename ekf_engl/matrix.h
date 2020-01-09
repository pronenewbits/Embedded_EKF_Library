/************************************************************************************
 *
 * Class Matrix
 *  Berisi kumpulan kode yang digunakan untuk melakukan operasi matrix.
 *
 *  Catatan:
 *    - Operasi indexing matrix dimulai dari 0, dengan format matrix[baris][kolom]
 *    - Data matrix disimpan dalam bentuk array 2 dimensi.
 *    - Setiap matrix yang menggunakan memori MATRIX_MAXIMUM_SIZE^2 (di variabel
 *       f32data), dan ukuran baris & kolom harus lebih kecil dari macro
 *       MATRIX_MAXIMUM_SIZE. Pendekatan ini digunakan untuk menghindari malloc
 *       (optimasi lebih lanjut bisa dilakukan untuk mengurangi penggunaan memori).
 * 
 * Class Matrix Versioning:
 *    v0.3_engl (2019-12-31), {PNb}:
 *      - Modifikasi dokumentasi kode buat orang asing.
 * 
 *    v0.3 (2019-12-25), {PNb}:
 *      - Menambahkan fungsi back subtitution untuk menyelesaikan permasalahan 
 *          persamaan linear Ax = B. Dengan A matrix segitiga atas & B vektor.
 *      - Memperbaiki bug pengecekan MATRIX_PAKAI_BOUND_CHECKING pada indexing kolom.
 *      - Menambahkan fungsi QR Decomposition (via Householder Transformation).
 *      - Menambahkan fungsi Householder Transformation.
 *      - Menghilangkan warning 'implicit conversion' untuk operasi pembandingan
 *          dengan float_prec_ZERO.
 *      - Menambahkan function overloading operasi InsertSubMatrix, untuk
 *          operasi insert dari SubMatrix ke SubMatrix.
 *      - Saat inisialisasi, matrix diisi nol (melalui vIsiHomogen(0.0)).
 *      - Menambahkan function overloading operator '/' dengan scalar.
 *
 *    v0.2 (2019-11-30), {PNb}:
 *      - Fungsi yang disupport:
 *          - Operator ==
 *          - Normalisasi matrix
 *          - Cholesky Decomposition
 *          - InsertSubMatrix
 *          - InsertVector
 *
 *    v0.1 (2019-11-29), {PNb}: 
 *      - Fungsi yang disupport:
 *          - Operasi matrix dasar
 *          - Invers
 *          - Cetak
 * 
 * See https://github.com/pronenewbits for more!
 *************************************************************************************/


#ifndef MATRIX_H
#define MATRIX_H

#include "konfig.h"

#if (SYSTEM_IMPLEMENTATION == SYSTEM_IMPLEMENTATION_PC)
    #include <iostream>
    #include <iomanip>      // std::setprecision

    using namespace std;
#elif (SYSTEM_IMPLEMENTATION == SYSTEM_IMPLEMENTATION_EMBEDDED_ARDUINO)
    #include <Wire.h>
#endif


class Matrix
{
public:
    Matrix(const int32_t _i32baris, const int32_t _i32kolom)
    {
        this->i32baris = _i32baris;
        this->i32kolom = _i32kolom;

        this->vIsiHomogen(0.0);
    }
    
    bool bCekMatrixValid() {
        /* Index terakhir untuk buffer jika ada kode yg buffer overflow 1 index */
        if ((this->i32baris > 0) && (this->i32baris < MATRIX_MAXIMUM_SIZE) && (this->i32kolom > 0) && (this->i32kolom < MATRIX_MAXIMUM_SIZE)) {
            return true;
        } else {
            return false;
        }
    }

    bool bCekMatrixPersegi() {
        return (this->i32baris == this->i32kolom);
    }
    
    int32_t i32getBaris() { return this->i32baris; }
    int32_t i32getKolom() { return this->i32kolom; }
    
    /* Ref: https://stackoverflow.com/questions/6969881/operator-overload */
    class Proxy {
    public:
        Proxy(float_prec* _array, int32_t _maxKolom) : _array(_array) { this->_maxKolom = _maxKolom; }

        /* Modifikasi agar lvalue modifiable, ref:
         * https://stackoverflow.com/questions/6969881/operator-overload#comment30831582_6969904
         * (I know this is so dirty, but it makes the code so FABULOUS :D)
         */
        float_prec & operator[](int32_t _kolom) {
            #if defined(MATRIX_PAKAI_BOUND_CHECKING)
                if (_kolom >= this->_maxKolom) {
                    #if (SYSTEM_IMPLEMENTATION == SYSTEM_IMPLEMENTATION_PC)
                        cout << "Matrix index out-of-bounds (kolom: " << _kolom << ")"<< endl;
                    #elif (SYSTEM_IMPLEMENTATION == SYSTEM_IMPLEMENTATION_EMBEDDED_ARDUINO)
                        Serial.println("Matrix index out-of-bounds kolom");
                    #else
                        /* Silent function */
                    #endif
                    while(1);
                }
            #endif
            return _array[_kolom];
        }
    private:
        float_prec* _array;
        int32_t _maxKolom;
    };
    Proxy operator[](int32_t _baris) {
        #if defined(MATRIX_PAKAI_BOUND_CHECKING)
            if (_baris >= this->i32baris) {
                #if (SYSTEM_IMPLEMENTATION == SYSTEM_IMPLEMENTATION_PC)
                    cout << "Matrix index out-of-bounds (baris: " << _baris << ")"<< endl;
                #elif (SYSTEM_IMPLEMENTATION == SYSTEM_IMPLEMENTATION_EMBEDDED_ARDUINO)
                    Serial.println("Matrix index out-of-bounds baris");
                #else
                    /* Silent function */
                #endif
                while(1);
            }
        #endif
        return Proxy(f32data[_baris], this->i32kolom);      /* Parsing data kolom untuk bound checking */
    }

    bool operator == (Matrix _pembanding) {
        if ((this->i32baris != _pembanding.i32baris) || (this->i32kolom != _pembanding.i32getKolom())) {
            return false;
        }

        for (int32_t _i = 0; _i < this->i32baris; _i++) {
            for (int32_t _j = 0; _j < this->i32kolom; _j++) {
                if (fabs((*this)[_i][_j] - _pembanding[_i][_j]) > float_prec(float_prec_ZERO)) {
                    return false;
                }
            }
        }
        return true;
    }

    Matrix operator + (Matrix _penjumlah) {
        Matrix _outp(this->i32baris, this->i32kolom);

        for (int32_t _i = 0; _i < this->i32baris; _i++) {
            for (int32_t _j = 0; _j < this->i32kolom; _j++) {
                _outp[_i][_j] = (*this)[_i][_j] + _penjumlah[_i][_j];
            }
        }
        return _outp;
    }

    Matrix operator - (Matrix _pengurang) {
        Matrix _outp(this->i32baris, this->i32kolom);

        for (int32_t _i = 0; _i < this->i32baris; _i++) {
            for (int32_t _j = 0; _j < this->i32kolom; _j++) {
                _outp[_i][_j] = (*this)[_i][_j] - _pengurang[_i][_j];
            }
        }
        return _outp;
    }

    Matrix operator * (Matrix _pengali) {
        Matrix _outp(this->i32baris, _pengali.i32kolom);

        for (int32_t _i = 0; _i < this->i32baris; _i++) {
            for (int32_t _j = 0; _j < _pengali.i32kolom; _j++) {
                _outp[_i][_j] = 0.0;
                for (int32_t _k = 0; _k < this->i32kolom; _k++) {
                    _outp[_i][_j] += ((*this)[_i][_k] * _pengali[_k][_j]);
                }
            }
        }
        return _outp;
    }

    Matrix operator * (float_prec _scalar) {
        Matrix _outp(this->i32baris, this->i32kolom);

        for (int32_t _i = 0; _i < this->i32baris; _i++) {
            for (int32_t _j = 0; _j < this->i32kolom; _j++) {
                _outp[_i][_j] = (*this)[_i][_j] * _scalar;
            }
        }
        return _outp;
    }

    Matrix operator / (float_prec _scalar) {
        Matrix _outp(this->i32baris, this->i32kolom);

        for (int32_t _i = 0; _i < this->i32baris; _i++) {
            for (int32_t _j = 0; _j < this->i32kolom; _j++) {
                _outp[_i][_j] = (*this)[_i][_j] / _scalar;
            }
        }
        return _outp;
    }

    void vRoundingElementToZero(const int32_t _i, const int32_t _j) {
        if (fabs((*this)[_i][_j]) < float_prec(float_prec_ZERO)) {
            (*this)[_i][_j] = 0.0;
        }
    }

    Matrix RoundingMatrixToZero() {
        for (int32_t _i = 0; _i < this->i32baris; _i++) {
            for (int32_t _j = 0; _j < this->i32kolom; _j++) {
                if (fabs((*this)[_i][_j]) < float_prec(float_prec_ZERO)) {
                    (*this)[_i][_j] = 0.0;
                }
            }
        }
        return (*this);
    }

    void vIsiHomogen(const float_prec _data) {
        for (int32_t _i = 0; _i < this->i32baris; _i++) {
            for (int32_t _j = 0; _j < this->i32kolom; _j++) {
                (*this)[_i][_j] = _data;
            }
        }
    }

    void vIsiNol() {
        this->vIsiHomogen(0.0);
    }

    void vIsiRandom(const int32_t _batasAtas, const int32_t _batasBawah) {
        for (int32_t _i = 0; _i < this->i32baris; _i++) {
            for (int32_t _j = 0; _j < this->i32kolom; _j++) {
                (*this)[_i][_j] = float_prec((rand() % (_batasAtas - _batasBawah + 1)) + _batasBawah);
            }
        }
    }

    void vIsiDiagonal(const float_prec _data) {
        for (int32_t _i = 0; _i < this->i32baris; _i++) {
            for (int32_t _j = 0; _j < this->i32kolom; _j++) {
                if (_i == _j) {
                    (*this)[_i][_j] = _data;
                } else {
                    (*this)[_i][_j] = 0.0;
                }
            }
        }
    }

    void vSetIdentitas() {
        this->vIsiDiagonal(1.0);
    }

    /* Masukkan vektor ke matrix pada posisi kolom _posKolom
     * Contoh: A = Matrix 3x3, B = Vector 3x1
     *
     *  C = A.InsertVector(B, 1);
     *
     *  A = [A00  A01  A02]     B = [B00]
     *      [A10  A11  A12]         [B10]
     *      [A20  A21  A22]         [B20]
     *
     *  C = [A00  B00  A02]
     *      [A10  B10  A12]
     *      [A20  B20  A22]
     */
    Matrix InsertVector(Matrix _Vector, const int32_t _posKolom) {
        Matrix _outp(this->i32kolom, this->i32baris);
        if ((_Vector.i32baris > this->i32baris) || (_Vector.i32kolom+_posKolom > this->i32kolom)) {
            /* Return false */
            _outp.i32baris = -1;
            _outp.i32kolom = -1;
            return _outp;
        }
        _outp = this->Salin();
        for (int32_t _i = 0; _i < _Vector.i32baris; _i++) {
            _outp[_i][_posKolom] = _Vector[_i][0];
        }
        return _outp;
    }

    /* Masukkan submatrix ke matrix pada posisi baris _posBaris & posisi kolom _posKolom
     * Contoh: A = Matrix 4x4, B = Matrix 2x3
     *
     *  C = A.InsertSubMatrix(B, 1, 1);
     *
     *  A = [A00  A01  A02  A03]    B = [B00  B01  B02]
     *      [A10  A11  A12  A13]        [B10  B11  B12]
     *      [A20  A21  A22  A23]
     *      [A30  A31  A32  A33]
     *
     *
     *  C = [A00  A01  A02  A03]
     *      [A10  B00  B01  B02]
     *      [A20  B10  B11  B12]
     *      [A30  A31  A32  A33]
     */
    Matrix InsertSubMatrix(Matrix _subMatrix, const int32_t _posBaris, const int32_t _posKolom) {
        Matrix _outp(this->i32kolom, this->i32baris);
        if (((_subMatrix.i32baris+_posBaris) > this->i32baris) || ((_subMatrix.i32kolom+_posKolom) > this->i32kolom)) {
            /* Return false */
            _outp.i32baris = -1;
            _outp.i32kolom = -1;
            return _outp;
        }
        _outp = this->Salin();
        for (int32_t _i = 0; _i < _subMatrix.i32baris; _i++) {
            for (int32_t _j = 0; _j < _subMatrix.i32kolom; _j++) {
                _outp[_i + _posBaris][_j + _posKolom] = _subMatrix[_i][_j];
            }
        }
        return _outp;
    }

    /* Masukkan submatrix dengan panjang _lenBaris & _lenKolom; ke submatrix pada posisi baris _posBaris & posisi kolom _posKolom
     * Contoh: A = Matrix 4x4, B = Matrix 2x3
     *
     *  C = A.InsertSubMatrix(B, 1, 1, 2, 2);
     *
     *  A = [A00  A01  A02  A03]    B = [B00  B01  B02]
     *      [A10  A11  A12  A13]        [B10  B11  B12]
     *      [A20  A21  A22  A23]
     *      [A30  A31  A32  A33]
     *
     *
     *  C = [A00  A01  A02  A03]
     *      [A10  B00  B01  A13]
     *      [A20  B10  B11  A23]
     *      [A30  A31  A32  A33]
     */
    Matrix InsertSubMatrix(Matrix _subMatrix, const int32_t _posBaris, const int32_t _posKolom, const int32_t _lenBaris, const int32_t _lenKolom) {
        Matrix _outp(this->i32kolom, this->i32baris);
        if (((_lenBaris+_posBaris) > this->i32baris) || ((_lenKolom+_posKolom) > this->i32kolom) || (_lenBaris > _subMatrix.i32baris) || (_lenKolom > _subMatrix.i32kolom)) {
            /* Return false */
            _outp.i32baris = -1;
            _outp.i32kolom = -1;
            return _outp;
        }
        _outp = this->Salin();
        for (int32_t _i = 0; _i < _lenBaris; _i++) {
            for (int32_t _j = 0; _j < _lenKolom; _j++) {
                _outp[_i + _posBaris][_j + _posKolom] = _subMatrix[_i][_j];
            }
        }
        return _outp;
    }

    /* Masukkan submatrix pada posisi baris _posBarisSub & posisi kolom _posKolomSub dengan panjang _lenBaris & _lenKolom;
     * ke submatrix pada posisi baris _posBaris & posisi kolom _posKolom.
     *
     * Contoh: A = Matrix 4x4, B = Matrix 2x3
     *
     *  C = A.InsertSubMatrix(B, 1, 1, 0, 1, 1, 2);
     *
     *  A = [A00  A01  A02  A03]    B = [B00  B01  B02]
     *      [A10  A11  A12  A13]        [B10  B11  B12]
     *      [A20  A21  A22  A23]
     *      [A30  A31  A32  A33]
     *
     *
     *  C = [A00  A01  A02  A03]
     *      [A10  B01  B02  A13]
     *      [A20  A21  A22  A23]
     *      [A30  A31  A32  A33]
     */
    Matrix InsertSubMatrix(Matrix _subMatrix, const int32_t _posBaris, const int32_t _posKolom,
                           const int32_t _posBarisSub, const int32_t _posKolomSub,
                           const int32_t _lenBaris, const int32_t _lenKolom) {
        Matrix _outp(this->i32kolom, this->i32baris);
        if (((_lenBaris+_posBaris) > this->i32baris) || ((_lenKolom+_posKolom) > this->i32kolom) ||
            ((_posBarisSub+_lenBaris) > _subMatrix.i32baris) || ((_posKolomSub+_lenKolom) > _subMatrix.i32kolom))
        {
            /* Return false */
            _outp.i32baris = -1;
            _outp.i32kolom = -1;
            return _outp;
        }
        _outp = this->Salin();
        for (int32_t _i = 0; _i < _lenBaris; _i++) {
            for (int32_t _j = 0; _j < _lenKolom; _j++) {
                _outp[_i + _posBaris][_j + _posKolom] = _subMatrix[_posBarisSub+_i][_posKolomSub+_j];
            }
        }
        return _outp;
    }

    Matrix Transpose() {
        Matrix _outp(this->i32kolom, this->i32baris);
        for (int32_t _i = 0; _i < this->i32baris; _i++) {
            for (int32_t _j = 0; _j < this->i32kolom; _j++) {
                _outp[_j][_i] = (*this)[_i][_j];
            }
        }
        return _outp;
    }

    bool bNormVector() {
        float_prec _normM = 0.0;
        for (int32_t _i = 0; _i < this->i32baris; _i++) {
            for (int32_t _j = 0; _j < this->i32kolom; _j++) {
                _normM = _normM + ((*this)[_i][_j] * (*this)[_i][_j]);
            }
        }
        if (_normM < float_prec(float_prec_ZERO)) {
            return false;
        }
        _normM = sqrt(_normM);
        for (int32_t _i = 0; _i < this->i32baris; _i++) {
            for (int32_t _j = 0; _j < this->i32kolom; _j++) {
                (*this)[_i][_j] /= _normM;
            }
        }
        return true;
    }
    
    Matrix Salin() {
        Matrix _outp(this->i32baris, this->i32kolom);
        for (int32_t _i = 0; _i < this->i32baris; _i++) {
            for (int32_t _j = 0; _j < this->i32kolom; _j++) {
                _outp[_i][_j] = (*this)[_i][_j];
            }
        }
        return _outp;
    }

    /* Melakukan operasi invers matrix dengan menggunakan algoritma Gauss-Jordan */
    Matrix Invers() {
        Matrix _outp(this->i32baris, this->i32kolom);
        Matrix _temp(this->i32baris, this->i32kolom);
        _outp.vSetIdentitas();
        _temp = this->Salin();


        /* Eliminasi Gauss... */
        for (int32_t _j = 0; _j < (_temp.i32baris)-1; _j++) {
            for (int32_t _i = _j+1; _i < _temp.i32baris; _i++) {
                if (fabs(_temp[_j][_j]) < float_prec(float_prec_ZERO)) {
                    // return false;    /* Matrix non-invertible */
                    _outp.i32baris = -1;
                    _outp.i32kolom = -1;
                    return _outp;
                }

                float_prec _tempfloat = _temp[_i][_j] / _temp[_j][_j];

                for (int32_t _k = 0; _k < _temp.i32kolom; _k++) {
                    _temp[_i][_k] -= (_temp[_j][_k] * _tempfloat);
                    _outp[_i][_k] -= (_outp[_j][_k] * _tempfloat);

                    _temp.vRoundingElementToZero(_i, _k);
                    _outp.vRoundingElementToZero(_i, _k);
                }

            }
        }

#if (1)
        /* Sampai sini seharusnya matrix _temp adalah matrix segitiga atas, tapi karena
         * keterbatasan kepresisian (rounding error), bisa jadi segitiga bawahnya
         * bukan 0 semua, jadikan 0 --> berguna untuk dekomposisi LU
         */
        for (int32_t _i = 1; _i < _temp.i32baris; _i++) {
            for (int32_t _j = 0; _j < _i; _j++) {
                _temp[_i][_j] = 0.0;
            }
        }
#endif


        /* Jordan... */
        for (int32_t _j = (_temp.i32baris)-1; _j > 0; _j--) {
            for (int32_t _i = _j-1; _i >= 0; _i--) {
                if (fabs(_temp[_j][_j]) < float_prec(float_prec_ZERO)) {
                    // return false;    /* Matrix non-invertible */
                    _outp.i32baris = -1;
                    _outp.i32kolom = -1;
                    return _outp;
                }

                float_prec _tempfloat = _temp[_i][_j] / _temp[_j][_j];
                _temp[_i][_j] -= (_temp[_j][_j] * _tempfloat);
                _temp.vRoundingElementToZero(_i, _j);

                for (int32_t _k = (_temp.i32baris - 1); _k >= 0; _k--) {
                    _outp[_i][_k] -= (_outp[_j][_k] * _tempfloat);
                    _outp.vRoundingElementToZero(_i, _k);
                }
            }
        }


        /* Normalisasi */
        for (int32_t _i = 0; _i < _temp.i32baris; _i++) {
            if (fabs(_temp[_i][_i]) < float_prec(float_prec_ZERO)) {
                // return false;    /* Matrix non-invertible */
                _outp.i32baris = -1;
                _outp.i32kolom = -1;
                return _outp;
            }

            float_prec _tempfloat = _temp[_i][_i];
            _temp[_i][_i] = 1.0;

            for (int32_t _j = 0; _j < _temp.i32baris; _j++) {
                _outp[_i][_j] /= _tempfloat;
            }
        }
        return _outp;
    }

    /* Melakukan operasi Dekomposisi Cholesky pada matrix dengan menggunakan algoritma Cholesky-Crout.
     *
     *      A = L*L'     ; A = matrix riil, positif definit, dan simetrik dengan ukuran MxM
     *
     *      L = A.CholeskyDec();
     *
     *      CATATAN! NOTE! Untuk menghemat komputansi, pengecekan matrix simetrik TIDAK dilakukan.
     *          Karena pemrosesan dilakukan pada segitiga kiri bawah dari matrix _A, maka
     *          diasumsikan segitiga atas dari _A juga merupakan simetrik dari segitiga bawah.
     *          (sebagai catatan, Scilab & MATLAB yang menggunakan Lapack routines DPOTRF
     *           memproses submatrix segitiga atas dari _A dan berperilaku kebalikan dengan
     *           fungsi ini, namun tetap valid secara matematis).
     *
     */
    Matrix CholeskyDec()
    {
        float_prec _tempFloat;

        Matrix _outp(this->i32baris, this->i32kolom);
        if (this->i32baris != this->i32kolom) {
            _outp.i32baris = -1;
            _outp.i32kolom = -1;
            return _outp;
        }
        _outp.vIsiHomogen(0.0);
        for (int32_t _j = 0; _j < this->i32kolom; _j++) {
            for (int32_t _i = _j; _i < this->i32baris; _i++) {
                _tempFloat = (*this)[_i][_j];
                if (_i == _j) {
                    for (int32_t _k = 0; _k < _j; _k++) {
                        _tempFloat = _tempFloat - (_outp[_i][_k] * _outp[_i][_k]);
                    }
                    if (_tempFloat < float_prec(float_prec_ZERO)) {
                        /* Matrix tidak positif definit */
                        _outp.i32baris = -1;
                        _outp.i32kolom = -1;
                        return _outp;
                    }
                    _outp[_i][_i] = sqrt(_tempFloat);
                } else {
                    for (int32_t _k = 0; _k < _j; _k++) {
                        _tempFloat = _tempFloat - (_outp[_i][_k] * _outp[_j][_k]);
                    }
                    if (fabs(_outp[_j][_j]) < float_prec(float_prec_ZERO)) {
                        /* Matrix tidak positif definit */
                        _outp.i32baris = -1;
                        _outp.i32kolom = -1;
                        return _outp;
                    }
                    _outp[_i][_j] = _tempFloat / _outp[_j][_j];
                }
            }
        }
        return _outp;
    }

    /* Melakukan operasi Householder Transformation pada matrix.
     *              out = HouseholderTransformQR(A, i, j)
     */
    Matrix HouseholderTransformQR(const int32_t _rowTransform, const int32_t _columnTransform)
    {
        float_prec _tempFloat;
        float_prec _xLen;
        float_prec _x1;
        float_prec _u1;
        float_prec _vLen2;

        Matrix _outp(this->i32baris, this->i32baris);
        Matrix _vectTemp(this->i32baris, 1);
        if ((_rowTransform >= this->i32baris) || (_columnTransform >= this->i32kolom)) {
            _outp.i32baris = -1;
            _outp.i32kolom = -1;
            return _outp;
        }

        /* sampai sini:
         *
         * _xLen    = ||x||            = sqrt(x1^2 + x2^2 + .. + xn^2)
         * _vLen2   = ||u||^2 - (u1^2) = x2^2 + .. + xn^2
         * _vectTemp= [0 0 0 .. x1=0 x2 x3 .. xn]'
         */
        _x1 = (*this)[_rowTransform][_columnTransform];
        _xLen = _x1*_x1;
        _vLen2 = 0.0;
        for (int32_t _i = _rowTransform+1; _i < this->i32baris; _i++) {
            _vectTemp[_i][0] = (*this)[_i][_columnTransform];

            _tempFloat = _vectTemp[_i][0] * _vectTemp[_i][0];
            _xLen  += _tempFloat;
            _vLen2 += _tempFloat;
        }
        _xLen = sqrt(_xLen);

        /* u1    = x1+(-sign(x1))*xLen */
        if (_x1 < 0.0) {
            _u1 = _x1+_xLen;
        } else {
            _u1 = _x1-_xLen;
        }


        /* Selesaikan vlen2 & tempHH */
        _vLen2 += (_u1*_u1);
        _vectTemp[_rowTransform][0] = _u1;

        if (fabs(_vLen2) < float_prec(float_prec_ZERO)) {
            /* vektor x sudah collinear dengan vektor basis e, kembalikan hasil = I */
            _outp.vSetIdentitas();
        } else {
            /* P = -2*(u1*u1')/v_len2 + I */
            /* PR TODO: Seharusnya bagian di bawah ini masih banyak yang bisa dioptimalkan */
            for (int32_t _i = 0; _i < this->i32baris; _i++) {
                _tempFloat = _vectTemp[_i][0];
                if (fabs(_tempFloat) > float_prec(float_prec_ZERO)) {
                    for (int32_t _j = 0; _j < this->i32baris; _j++) {
                        if (fabs(_vectTemp[_j][0]) > float_prec(float_prec_ZERO)) {
                            _outp[_i][_j] = _vectTemp[_j][0];
                            _outp[_i][_j] = _outp[_i][_j] * _tempFloat;
                            _outp[_i][_j] = _outp[_i][_j] * (-2.0/_vLen2);
                        }
                    }
                }
                _outp[_i][_i] = _outp[_i][_i] + 1.0;
            }
        }
        return _outp;
    }

    /* Melakukan operasi QR decomposition pada matrix.
     *      Melakukan operasi Dekomposisi QR pada matrix A dengan menggunakan fungsi
     *          MatrixOp_HouseholderTransform.
     *                      A = Q*R
     */
    bool QRDec(Matrix &Q, Matrix &R)
    {
        Matrix Qn(Q.i32baris, Q.i32kolom);
        if ((this->i32baris < this->i32kolom) || (!Q.bCekMatrixPersegi()) || (Q.i32baris != this->i32baris) || (R.i32baris != this->i32baris) || (R.i32kolom != this->i32kolom)) {
            Q.i32baris = -1;
            Q.i32kolom = -1;
            R.i32baris = -1;
            R.i32kolom = -1;
            return false;
        }
        R = (*this);
        Q.vSetIdentitas();
        for (int32_t _i = 0; (_i < (this->i32baris - 1)) && (_i < this->i32kolom-1); _i++) {
            Qn  = R.HouseholderTransformQR(_i, _i);
            Q   = Qn * Q;
            R   = Qn * R;
        }
        Q.RoundingMatrixToZero();
        R.RoundingMatrixToZero();
        Q = Q.Transpose();
        return true;
    }


//    /* Melakukan operasi Forward-subtitution pada matrix triangular A & matrix kolom B.
//     *                      Ax = B
//     *
//     *  Untuk menghemat komputansi, matrix A tidak dilakukan pengecekan triangular
//     * (diasumsikan sudah lower-triangular).
//     */
//    Matrix ForwardSubtitution(Matrix &A, Matrix &B)
//    {
//        Matrix _outp(A.i32baris, 1);
//        if ((A.i32baris != A.i32kolom) || (A.i32baris != B.i32baris)) {
//            _outp.i32baris = -1;
//            _outp.i32kolom = -1;
//            return _outp;
//        }

//        for (int32_t _i = 0; _i < A.i32baris; _i++) {
//            _outp[_i][0] = B[_i][0];
//            for (int32_t _j = 0; _j < _i; _j++) {
//                _outp[_i][0] = _outp[_i][0] - A[_i][_j]*_outp[_j][0];
//            }
//            if (fabs(A[_i][_i]) < float_prec(float_prec_ZERO)) {
//                _outp.i32baris = -1;
//                _outp.i32kolom = -1;
//                return _outp;
//            }
//            _outp[_i][0] = _outp[_i][0] / A[_i][_i];
//        }
//        return _outp;
//    }


    /* Melakukan operasi back-subtitution pada matrix triangular A & matrix kolom B.
     *                      Ax = B
     *
     *  Untuk menghemat komputansi, matrix A tidak dilakukan pengecekan triangular
     * (diasumsikan sudah upper-triangular).
     */
    Matrix BackSubtitution(Matrix &A, Matrix &B)
    {
        Matrix _outp(A.i32baris, 1);
        if ((A.i32baris != A.i32kolom) || (A.i32baris != B.i32baris)) {
            _outp.i32baris = -1;
            _outp.i32kolom = -1;
            return _outp;
        }

        for (int32_t _i = A.i32kolom-1; _i >= 0; _i--) {
            _outp[_i][0] = B[_i][0];
            for (int32_t _j = _i + 1; _j < A.i32kolom; _j++) {
                _outp[_i][0] = _outp[_i][0] - A[_i][_j]*_outp[_j][0];
            }
            if (fabs(A[_i][_i]) < float_prec(float_prec_ZERO)) {
                _outp.i32baris = -1;
                _outp.i32kolom = -1;
                return _outp;
            }
            _outp[_i][0] = _outp[_i][0] / A[_i][_i];
        }

        return _outp;
    }

#if (SYSTEM_IMPLEMENTATION == SYSTEM_IMPLEMENTATION_PC)
    void vCetak() {
        for (int32_t _i = 0; _i < this->i32baris; _i++) {
            cout << "[ ";
            for (int32_t _j = 0; _j < this->i32kolom; _j++) {
                cout << std::fixed << std::setprecision(3) << (*this)[_i][_j] << " ";
            }
            cout << "]";
            cout << endl;
        }
        cout << endl;
    }
    void vCetakFull() {
        for (int32_t _i = 0; _i < this->i32baris; _i++) {
            cout << "[ ";
            for (int32_t _j = 0; _j < this->i32kolom; _j++) {
                cout << resetiosflags( ios::fixed | ios::showpoint ) << (*this)[_i][_j] << " ";
            }
            cout << "]";
            cout << endl;
        }
        cout << endl;
    }
#elif (SYSTEM_IMPLEMENTATION == SYSTEM_IMPLEMENTATION_EMBEDDED_ARDUINO)
    void vCetak() {
        char _bufSer[10];
        for (int32_t _i = 0; _i < this->i32baris; _i++) {
            Serial.print("[ ");
            for (int32_t _j = 0; _j < this->i32kolom; _j++) {
                snprintf(_bufSer, sizeof(_bufSer)-1, "%2.2f ", (*this)[_i][_j]);
                Serial.print(_bufSer);
            }
            Serial.println("]");
        }
        Serial.println("");
    }
#else
    #warning("Fungsi Matrix.vCetak() tidak berfungsi");
    
    void vCetak() {}     /* Silent function */
#endif

private:
    int32_t i32baris;
    int32_t i32kolom;
    float_prec f32data[MATRIX_MAXIMUM_SIZE][MATRIX_MAXIMUM_SIZE] = {{0}};
};

#endif // MATRIX_H
