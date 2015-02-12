/*! \file scfloop.h
    \brief SCFを行うクラスの宣言

    Copyright ©  2015 @dc1394 All Rights Reserved.
*/

#ifndef _SCFLOOP_H_
#define _SCFLOOP_H_

#pragma once

#include "diffdata.h"

namespace schrac {
	class ScfLoop final {
        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            唯一のコンストラクタ
            \param arg インプットファイル名とTBBを使用するかどうかのstd::pair
        */
        explicit ScfLoop(std::pair<std::string, bool> const & arg);
        
        //! A destructor.
        /*!
            何もしないデストラクタ
        */
        ~ScfLoop()
        {
        }

        // #region メンバ関数
        
        bool check_converge();

        // #endregion メンバ関数

        // #region メンバ変数

    private:
        //!  A private member variable (constant).
        /*!
            微分方程式のデータオブジェクト
        */
        std::shared_ptr<DiffData> const pdiffdata_;

        // #endregion メンバ変数
	};

	/*inline void SCFLoop::setExp() const
	{
		std::cout.setf(std::ios::fixed, std::ios::floatfield);
		std::cout << std::setprecision(
			boost::numeric_cast<const std::streamsize>(std::fabs(std::log10(pdata_->eps))) - 2);
	}*/
}

#endif  // _SCFLOOP_H_
