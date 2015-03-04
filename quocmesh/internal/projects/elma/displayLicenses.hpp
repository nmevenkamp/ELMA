#ifndef __DISPLAYLICENSES_H
#define __DISPLAYLICENSES_H

#include <QApplication>
#include <QMenu>
#include <QDialog>
#include <QMessageBox>
#include <QSpacerItem>
#include <QGridLayout>

#include <sstream>
#ifdef USE_LIB_PNG
#include <png.h>
#endif

class QAboutExternalLibrariesMenu : public QMenu {
  Q_OBJECT
public:
  QAction *actionAboutBLAS;
  QAction *actionAboutBoost;
  QAction *actionAboutKissFFT;
  QAction *actionAboutLAPACK;
  QAction *actionAboutLibpng;
  QAction *actionAboutLibTIFF;
  QAction *actionAboutWavelet1d;
  QAction *actionAboutZlib;
  
  QAboutExternalLibrariesMenu ( QWidget * parent = 0 ) : QMenu ( parent ) {
    /*
     * BEGIN: About external libraries
     */    
    // About BLAS
    actionAboutBLAS = new QAction(this);
    actionAboutBLAS->setObjectName(QString::fromUtf8("actionAboutBLAS"));
    actionAboutBLAS->setText(QApplication::translate("QuocImageViewer", "About BLAS", 0));
    this->addAction(actionAboutBLAS);
    QObject::connect(actionAboutBLAS, SIGNAL(triggered()), this, SLOT(on_actionAboutBLAS_triggered()));
    
    // About Boost
    actionAboutBoost = new QAction(this);
    actionAboutBoost->setObjectName(QString::fromUtf8("actionAboutBOOST"));
    actionAboutBoost->setText(QApplication::translate("QuocImageViewer", "About Boost", 0));
    this->addAction(actionAboutBoost);
    QObject::connect(actionAboutBoost, SIGNAL(triggered()), this, SLOT(on_actionAboutBoost_triggered()));
    
    // About LAPACK
    actionAboutLAPACK = new QAction(this);
    actionAboutLAPACK->setObjectName(QString::fromUtf8("actionAboutLAPACK"));
    actionAboutLAPACK->setText(QApplication::translate("QuocImageViewer", "About LAPACK", 0));
    this->addAction(actionAboutLAPACK);
    QObject::connect(actionAboutLAPACK, SIGNAL(triggered()), this, SLOT(on_actionAboutLAPACK_triggered()));
    
    // About libpng
    actionAboutLibpng = new QAction(this);
    actionAboutLibpng->setObjectName(QString::fromUtf8("actionAboutLibpng"));
    actionAboutLibpng->setText(QApplication::translate("QuocImageViewer", "About libpng", 0));
    this->addAction(actionAboutLibpng);
    QObject::connect(actionAboutLibpng, SIGNAL(triggered()), this, SLOT(on_actionAboutLibpng_triggered()));
    
    // About LibTIFF
    actionAboutLibTIFF = new QAction(this);
    actionAboutLibTIFF->setObjectName(QString::fromUtf8("actionAboutLibTIFF"));
    actionAboutLibTIFF->setText(QApplication::translate("QuocImageViewer", "About LibTIFF", 0));
    this->addAction(actionAboutLibTIFF);
    QObject::connect(actionAboutLibTIFF, SIGNAL(triggered()), this, SLOT(on_actionAboutLibTIFF_triggered()));
    
    // About Kiss FFT
    actionAboutKissFFT = new QAction(this);
    actionAboutKissFFT->setObjectName(QString::fromUtf8("actionAboutKissFFT"));
    actionAboutKissFFT->setText(QApplication::translate("QuocImageViewer", "About Kiss FFT", 0));
    this->addAction(actionAboutKissFFT);
    QObject::connect(actionAboutKissFFT, SIGNAL(triggered()), this, SLOT(on_actionAboutKissFFT_triggered()));
    
    // About wavelet1d
    actionAboutWavelet1d = new QAction(this);
    actionAboutWavelet1d->setObjectName(QString::fromUtf8("actionAboutWavelet1d"));
    actionAboutWavelet1d->setText(QApplication::translate("QuocImageViewer", "About wavelet1d", 0));
    this->addAction(actionAboutWavelet1d);
    QObject::connect(actionAboutWavelet1d, SIGNAL(triggered()), this, SLOT(on_actionAboutWavelet1d_triggered()));
    
    // About zlib
    actionAboutZlib = new QAction(this);
    actionAboutZlib->setObjectName(QString::fromUtf8("actionAboutZlib"));
    actionAboutZlib->setText(QApplication::translate("QuocImageViewer", "About zlib", 0));
    this->addAction(actionAboutZlib);
    QObject::connect(actionAboutZlib, SIGNAL(triggered()), this, SLOT(on_actionAboutZlib_triggered()));
    /*
     * END: About external libraries
     */
  }
  
protected slots:
  void on_actionAboutBLAS_triggered () {
    QMessageBox msgBox(this);
    msgBox.setWindowTitle("About BLAS (Basic Linear Algebra Subprograms)");
    msgBox.setTextFormat(Qt::RichText);
    msgBox.setText("<font size='4'><b>About BLAS</b><br><br></font><font size='2'>Basic Linear Algebra Subprograms<br><br>The reference BLAS is a freely-available software package. It is available from netlib via anonymous ftp and the World Wide Web. Thus, it can be included in commercial software packages (and has been). We only ask that proper credit be given to the authors.<br><br>Like all software, it is copyrighted. It is not trademarked, but we do ask the following:<ul><li>If you modify the source for these routines we ask that you change the name of the routine and comment the changes made to the original.</li><li>We will gladly answer any questions regarding the software. If a modification is done, however, it is the responsibility of the person who modified the routine to provide support.</li></ul><br>For further information, we refer to: <a href='http://www.netlib.org/blas/'>http://www.netlib.org/blas/</a></font>");
    QSpacerItem* horizontalSpacer = new QSpacerItem(400, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);
    QGridLayout* layout = (QGridLayout*)msgBox.layout();
    layout->addItem(horizontalSpacer, layout->rowCount(), 0, 1, layout->columnCount());
    msgBox.exec();
  }
  
  void on_actionAboutBoost_triggered () {
    QMessageBox msgBox(this);
    msgBox.setWindowTitle("About Boost");
    msgBox.setTextFormat(Qt::RichText);
    msgBox.setText("<font size='4'><b>About Boost</b><br><br></font><font size='2'>Boost Software License - Version 1.0 - August 17th, 2003<br><br>Permission is hereby granted, free of charge, to any person or organization obtaining a copy of the software and accompanying documentation covered by this license (the \"Software\") to use, reproduce, display, distribute, execute, and transmit the Software, and to prepare derivative works of the Software, and to permit third-parties to whom the Software is furnished to do so, all subject to the following:<br><br>The copyright notices in the Software and this entire statement, including the above license grant, this restriction and the following disclaimer, must be included in all copies of the Software, in whole or in part, and all derivative works of the Software, unless such copies or derivative works are solely in the form of machine-executable object code generated by a source language processor.<br><br> THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.<br><br>For further information, we refer to: <a href='http://www.boost.org'>http://www.boost.org</a></font>");
    QSpacerItem* horizontalSpacer = new QSpacerItem(400, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);
    QGridLayout* layout = (QGridLayout*)msgBox.layout();
    layout->addItem(horizontalSpacer, layout->rowCount(), 0, 1, layout->columnCount());
    msgBox.exec();
  }
  
  void on_actionAboutLAPACK_triggered () {
    QMessageBox msgBox(this);
    msgBox.setWindowTitle("About LAPACK (Linear Algebra PACKage)");
    msgBox.setTextFormat(Qt::RichText);
    msgBox.setText("<font size='4'><b>About LAPACK</b><br><br></font><font size='2'>Linear Algebra PACKage<br><table><tr><td>Copyright (c) 1992-2013</td><td>The University of Tennessee and The University of Tennessee Research Foundation. All rights reserved.</td></tr><tr><td>Copyright (c) 2000-2013</td><td>The University of California Berkeley. All rights reserved.</td></tr><tr><td>Copyright (c) 2006-2013</td><td>The University of Colorado Denver. All rights reserved.</td></tr></table><br><br>$COPYRIGHT$<br><br>Additional copyrights may follow<br><br>$HEADER$<br><br>Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:<ul><li>Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.</li><li>Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer listed in this license in the documentation and/or other materials provided with the distribution.</li><li>Neither the name of the copyright holders nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.</li></ul><br>The copyright holders provide no reassurances that the source code provided does not infringe any patent, copyright, or any other intellectual property rights of third parties.  The copyright holders disclaim any liability to any recipient for claims brought against recipient by any third party for infringement of that parties intellectual property rights.<br><br>THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.<br><br>For further information, we refer to: <a href='http://www.netlib.org/lapack/'>http://www.netlib.org/lapack/</a></font>");
    QSpacerItem* horizontalSpacer = new QSpacerItem(600, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);
    QGridLayout* layout = (QGridLayout*)msgBox.layout();
    layout->addItem(horizontalSpacer, layout->rowCount(), 0, 1, layout->columnCount());
    msgBox.exec();
  }
  
  void on_actionAboutLibpng_triggered () {
    QMessageBox msgBox(this);
    msgBox.setWindowTitle("About libpng");
    msgBox.setTextFormat(Qt::RichText);
    std::stringstream ss;
    ss << "<font size='4'><b>About libpng</b><br><br></font><font size='2'>";
#ifdef USE_LIB_PNG
    ss << png_get_copyright(NULL);
#else
    ss << "libpng was not used for compiling this executable";
#endif
    ss << "<br><br>For further information, we refer to: <a href='http://www.libpng.org/pub/png/libpng.html'>http://www.libpng.org/pub/png/libpng.html</a></font>";
    msgBox.setText( ss.str ( ).c_str ( ) );
    QSpacerItem* horizontalSpacer = new QSpacerItem(400, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);
    QGridLayout* layout = (QGridLayout*)msgBox.layout();
    layout->addItem(horizontalSpacer, layout->rowCount(), 0, 1, layout->columnCount());
    msgBox.exec();
  }
  
  void on_actionAboutLibTIFF_triggered () {
    QMessageBox msgBox(this);
    msgBox.setWindowTitle("About LibTIFF");
    msgBox.setTextFormat(Qt::RichText);
    msgBox.setText("<font size='4'><b>About LibTIFF</b><br><br></font><font size='2'>Copyright (c) 1988-1997 Sam Leffler<br>Copyright (c) 1991-1997 Silicon Graphics, Inc.<br><br>Permission to use, copy, modify, distribute, and sell this software and its documentation for any purpose is hereby granted without fee, provided that (i) the above copyright notices and this permission notice appear in all copies of the software and related documentation, and (ii) the names of Sam Leffler and Silicon Graphics may not be used in any advertising or publicity relating to the software without the specific, prior written permission of Sam Leffler and Silicon Graphics.<br><br>THE SOFTWARE IS PROVIDED \"AS-IS\" AND WITHOUT WARRANTY OF ANY KIND, EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.<br><br>IN NO EVENT SHALL SAM LEFFLER OR SILICON GRAPHICS BE LIABLE FOR ANY SPECIAL, INCIDENTAL, INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY KIND, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER OR NOT ADVISED OF THE POSSIBILITY OF DAMAGE, AND ON ANY THEORY OF LIABILITY, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.<br><br>For further information, we refer to: <a href='http://www.libtiff.org'>http://www.libtiff.org</a></font>");
    QSpacerItem* horizontalSpacer = new QSpacerItem(400, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);
    QGridLayout* layout = (QGridLayout*)msgBox.layout();
    layout->addItem(horizontalSpacer, layout->rowCount(), 0, 1, layout->columnCount());
    msgBox.exec();
  }
  
  void on_actionAboutKissFFT_triggered () {
    QMessageBox msgBox(this);
    msgBox.setWindowTitle("About Kiss FFT");
    msgBox.setTextFormat(Qt::RichText);
    msgBox.setText("<font size='4'><b>About Kiss FFT</b><br><br></font><font size='2'>Copyright (c) 2003-2010 Mark Borgerding<br><br>All rights reserved.<br><br>Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:<ul><li>Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.</li><li>Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.</li><li>Neither the author nor the names of any contributors may be used to endorse or promote products derived from this software without specific prior written permission.</li></ul><br>THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A ARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.<br><br>For further information, we refer to: <a href='http://sourceforge.net/projects/kissfft/'>http://sourceforge.net/projects/kissfft/</a></font>");
    QSpacerItem* horizontalSpacer = new QSpacerItem(400, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);
    QGridLayout* layout = (QGridLayout*)msgBox.layout();
    layout->addItem(horizontalSpacer, layout->rowCount(), 0, 1, layout->columnCount());
    msgBox.exec();
  }
  
  void on_actionAboutWavelet1d_triggered () {
    QMessageBox msgBox(this);
    msgBox.setWindowTitle("About wavelet1d");
    msgBox.setTextFormat(Qt::RichText);
    msgBox.setText("<font size='4'><b>About wavelet1d</b><br><br></font><font size='2'>C++ 1D/2D DWT Implementation for Win32 and Linux<br><br>Copyright (C) &lt;2011&gt; by &lt;Rafat Hussain&gt;<br><br>Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the \"Software\"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:<br><br>The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.<br><br> THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.<br><br>For further information, we refer to: <a href='http://code.google.com/p/wavelet1d/'>http://code.google.com/p/wavelet1d/</a></font>");
    msgBox.exec();
  }
  
  void on_actionAboutZlib_triggered () {
    QMessageBox msgBox(this);
    msgBox.setWindowTitle("About zlib");
    msgBox.setTextFormat(Qt::RichText);
    msgBox.setText("<font size='4'><b>About zlib</b><br><br></font><font size='2'>A Massively Spiffy Yet Delicately Unobtrusive Compression Library<br><br>Copyright (C) 1995-2013 Jean-loup Gailly and Mark Adler<br><br>This software is provided 'as-is', without any express or implied warranty.  In no event will the authors be held liable for any damages arising from the use of this software.<br><br> Permission is granted to anyone to use this software for any purpose, including commercial applications, and to alter it and redistribute it freely, subject to the following restrictions:<ol><li>The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.</li><li>Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.</li><li>This notice may not be removed or altered from any source distribution.</li><br><table><tr><td>Jean-loup Gailly</td><td>Mark Adler</td></tr><tr><td>jloup@gzip.org</td><td>madler@alumni.caltech.edu</td></tr></table><br><br>For further information, we refer to: <a href='http://www.zlib.net'>http://www.zlib.net</a></font>");
    QSpacerItem* horizontalSpacer = new QSpacerItem(400, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);
    QGridLayout* layout = (QGridLayout*)msgBox.layout();
    layout->addItem(horizontalSpacer, layout->rowCount(), 0, 1, layout->columnCount());
    msgBox.exec();
  }
};

#endif // __DISPLAYLICENSES_H
