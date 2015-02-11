/*! \file checkpoint.h
    \brief 時間計測のためのクラスの宣言

    Copyright ©  2014 @dc1394 All Rights Reserved.
*/
#ifndef _CHECKPOINT_H_
#define _CHECKPOINT_H_

#pragma once

#include <cstdint>              // for std::int32_t, std::int64_t
#include <memory>               // for std::unique_ptr
#include <utility>              // for std::pair

namespace checkpoint {
    //! A class.
    /*!
        時間計測のためのクラス
    */
	class CheckPoint final {
        // #region クラスの前方宣言

        struct Timestamp;
        struct CheckPointFastImpl;

        // #endregion クラスの前方宣言

        // #region クラス内クラスの宣言と実装

        template <typename T>
        struct fastpimpl_deleter {
            void operator()(T * p) const {
                FastArenaObject<sizeof(CheckPointFastImpl)>::
                    operator delete(reinterpret_cast<void *>(p));
            }
        };

        // #endregion クラス内クラスの宣言と実装

    public:
        // #region コンストラクタ・デストラクタ

        //! A constructor.
        /*!
            デフォルトコンストラクタかつ唯一のコンストラクタ
        */
        CheckPoint();

        //! A destructor.
        /*!
        */
        ~CheckPoint();

        // #endregion コンストラクタ・デストラクタ

        // #region メンバ関数

        //! A public member function.
        /*!
            チェックポイントを設定する
            \param line 行数
            \param action チェックポイントの名称
        */
        void checkpoint(char const * action, std::int32_t line);

        //! A public member function.
        /*!
            直前のチェックポイントから計測した、経過時間を表示する
        */
        void checkpoint_print() const;

        //! A public member function.
        /*!
            最初のチェックポイントから最後のチェックポイント
            までの経過時間を返す
            \return 経過時間
        */
        double totalpassageoftime() const;

        // #endregion メンバ関数 

    private:
        // #region メンバ変数

        //! A private member variable (constant).
        /*!
            チェックポイントの状態へのスマートポインタ
        */
        const std::unique_ptr<CheckPointFastImpl, fastpimpl_deleter<CheckPointFastImpl>> cfp;

        // #endregion メンバ変数

        // #region 禁止されたコンストラクタ・メンバ関数

        //! A private copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
        */
		CheckPoint(CheckPoint const &) = delete;

        //! operator=() (deleted).
        /*!
            operator=()の宣言（禁止）
            \param コピー元のオブジェクト
            \return コピー元のオブジェクト
        */
		CheckPoint & operator=(CheckPoint const &) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
	};

#ifdef _WIN32
    void usedmem();
#endif
}

#endif  // _CHECKPOINT_H_
