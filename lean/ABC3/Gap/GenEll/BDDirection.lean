import ABC3.Found.GenEll.BDClass

/-!
# Gap —— [GenEll] `≲` の向きが、定義と用法で食い違う

★**これは「読めなかった」でも「証明できなかった」でもない。**
★**印字された定義と、同じ論文の中の用法が、互いに逆を向いている**ことを特定した段階である。

## ★事実 1 —— 定義(物理 p.5、260 dpi 目視確認 2026-08-16)

原文 (GenEll p.5):
> if there exists a ["constant"] C ∈ R such that β(x) − α(x) ≤ C (respectively,

すなわち `α ≲_F β` は **`β(x) − α(x) ≤ C`**、言い換えれば **`β ≤ α + C`**。
★**左辺が上から押さえるのではなく、右辺が押さえられる**。

## ★事実 2 —— 用法(物理 p.9・p.11・p.6)

| 箇所 | 印字 | 表題・文脈が求める向き |
|---|---|---|
| `Proposition 1.6` | `log-cond_D ≲ ht_D` | 表題 "Conductor **Bounded by** the Height" ⇒ `log-cond ≤ ht + C` |
| `Theorem 2.1`, (i)(ii) | `ht ≲ (1+ϵ)(log-diff + log-cond)` | abc 予想 ⇒ `ht ≤ (1+ϵ)(…) + C` |
| `Proposition 1.4`, (ii) | `ht_L̄ ≳ 0` | 「大域切断で生成されるなら高さは**下に**有界」 |

★**3 箇所すべてで、定義どおりに読むと向きが逆になる。**

## ★事実 3 —— 向きの違いは無害ではない(機械的な証拠)

`Found/GenEll/BDClass.lean` の `bdle_ne_bdge` が
**`BDle α β ∧ ¬ BDge α β` なる `α, β` を実際に構成している**
(`α = 0`、`β(n) = −n`)。
★すなわち **`≲` と `≳` は同じ関係ではない**。取り違えれば主張は別物になる。

## ★★分類は ① `modelError` である —— ③ を名乗らない

★既定は `modelError`(①)であり、`sourceGap`(③)は最後の手段である。
本件を ① とする理由:

1. ★**この食い違いは「印刷上の慣習」で説明が付く可能性が残っている。**
   `≲` を「左辺の**誤差**が右辺で押さえられる」と読む流儀
   (`α ≲ β` ⇔ 「α は β より本質的に小さくない」)は文献に存在する。
2. ★**我々は他の Mochizuki 論文での同記法の用例をまだ数えていない。**
   `[IUTchIV]` は同じ `≲` を使う。**そちらの用法を突き合わせるまでは、
   こちらの読み違いである可能性を排除できない。**

★**`falsifier` を書けないうちは ③ を名乗ってはならない**——本ファイルはその規律に従う。

## ★我々が取った措置

★**印字どおりに写した。**「表題が求める向き」に直して写すことはしていない。
- `Found/GenEll/BDClass.lean` —— `BDle` は逐語(`β(x) − α(x) ≤ C`)
- `Skeleton/GenEll/Section1.lean` —— `Proposition 1.4, (ii)` / `Proposition 1.6`
- `Skeleton/GenEll/Section2.lean` —— `Theorem 2.1`, (i)(ii)

★**直してしまうと、この食い違いは二度と検出できなくなる。**
写しを原文に合わせ、食い違いは記録に置く——この分業が本ファイルの目的である。
-/

namespace ABC3.Gap.GenEll

open ABC3.Meta ABC3.Found.GenEll

/-- **[GenEll] Definition 1.2, (ii)** の `≲` の向きについての記録。 -/
def bdDirection : GapRecord :=
  { source :=
      { paper := "GenEll", pdfPage := 5, item := "Definition 1.2, (ii)",
        sectionId := "genell-def-1-2-ii" }
    classification := .modelError
    falsifier :=
      "★① から動かす条件は 2 つ。(a) [IUTchIV] を含む Mochizuki の他論文で同じ `≲` の用例を数え、そこでも定義と用法が逆であることを示せば、『こちらの読み違い』の線が消えて ③(sourceGap)へ上がる。(b) 逆に、`α ≲ β` を『α は β より本質的に小さくない』と読む流儀の下では p.5 の定義と p.6・p.9・p.11 の用法が両立することを示せれば、本記録は ①(我々のモデル化の誤り)として閉じる。★どちらもまだ測っていない。" }

/-- ★向きの違いが**無害でない**ことの機械的な証拠。

`Found/GenEll/BDClass.lean` の `bdle_ne_bdge` を、この Gap の witness としてここに束ねる。
★**記録だけでは主張にならない**——反例が要る。 -/
theorem bdDirection_matters :
    ∃ α β : ℕ → ℝ, BDle α β ∧ ¬ BDge α β :=
  bdle_ne_bdge

end ABC3.Gap.GenEll
