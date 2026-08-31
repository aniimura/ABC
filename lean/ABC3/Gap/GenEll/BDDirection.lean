import ABC3.Found.GenEll.BDClass
import ABC3.Check.GenEll.Prop17Direction

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

## ★★★★★★★★★★事実 2b —— 4 例目は「逆」ではなく**偽**である（2026-08-29 追記）

| 箇所 | 印字 | 印字どおりに読むと |
|---|---|---|
| `Proposition 1.7`, (i) 右 | `log-diff の差 ≲ (1−1/e)·log-cond_E` | ★★**偽になる** |

★`Check/GenEll/Prop17Direction.lean` の `prop_1_7_printed_direction_false` で
**機械検証した**。反例は `φ = id`（`K = F`、`D ≔ E_red`）である:

* 条件 (b)(c)(d) はすべて満たされる（恒等射、`e_P = 1` は任意の `e` を割る）
* `log-diff(K) − log-diff(F) = 0`（同じ体）
* `log-cond_D = log-cond_E`（`log-cond` は根基で定義されている——`§9-966`）
* ★我々が証明した**左の等式**（`log-cond_E − log-cond_D = log-diff の差`）は成り立つ
* ★★ところが印字どおりの右の `≲` は `(1−1/e)·log-cond_E(x) ≤ C` を要求する
  ——導手は点にわたって非有界だから、そのような `C` は無い

★★★★★**通常の読み**（`log-diff の差 ≤ (1−1/e)·log-cond_E + C`）なら成り立つ
（`prop_1_7_ordinary_direction_holds`）。
★★これは「表題・文脈と逆」ではなく、**主張そのものが破れる**初めての例である。

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

## ★★★★★2026-08-29 —— `falsifier` は書けた。それでも ③ を名乗らない

事実 2b（`Proposition 1.7`, (i) の右の `≲` が印字どおりでは**偽**）は、
「印字どおりに読むと論文の主張が破れる」ことを**機械で**示している。
★これは本ファイルが待っていた `falsifier` の**片方**である。

★★それでも分類を ③ に上げないのは、上の理由 1・2 が**そのまま生きている**からである
——反例が示すのは「印字どおりの読みでは偽」であって、
「原典に穴がある」ことではない。★誤植の可能性と、`≲` の流儀を我々が読み違えている
可能性の**両方**が残る。★★★③ へ上げる条件は上の `falsifier` フィールドの (a)
（他論文での用例を数えること）のままである。

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
      "★① から動かす条件は 2 つ。(a) [IUTchIV] を含む Mochizuki の他論文で同じ `≲` の用例を数え、そこでも定義と用法が逆であることを示せば、『こちらの読み違い』の線が消えて ③(sourceGap)へ上がる。(b) 逆に、`α ≲ β` を『α は β より本質的に小さくない』と読む流儀の下では p.5 の定義と p.6・p.9・p.11 の用法が両立することを示せれば、本記録は ①(我々のモデル化の誤り)として閉じる。★★★2026-08-29: (a) の**前半**が取れた——`Check/GenEll/Prop17Direction.lean` の `prop_1_7_printed_direction_false` が『印字どおりに読むと Proposition 1.7, (i) の右の ≲ は偽』を機械検証した(φ = id の反例)。★用例表の 4 例目にして、初めて『文脈と逆』ではなく『偽になる』例である。★★それでも ③ へは上げない——残るのは他論文での用例を数える段であり、誤植の可能性と読み違いの可能性が両方残っているからである。" }

/-- ★向きの違いが**無害でない**ことの機械的な証拠。

`Found/GenEll/BDClass.lean` の `bdle_ne_bdge` を、この Gap の witness としてここに束ねる。
★**記録だけでは主張にならない**——反例が要る。 -/
theorem bdDirection_matters :
    ∃ α β : ℕ → ℝ, BDle α β ∧ ¬ BDge α β :=
  bdle_ne_bdge

/-- ★★★★★★★★★★**向きの食い違いが「主張を破る」ことの機械的な証拠**（2026-08-29）。

`Check/GenEll/Prop17Direction.lean` の反例を、この Gap の witness としてここに束ねる。
★`bdDirection_matters` は「`≲` と `≳` は別の関係である」ことしか言っていない。
★★**本定理は『印字どおりに読むと [GenEll] Proposition 1.7, (i) が偽になる』と言う**
——用例表の 4 例目にして、初めて「文脈と逆」ではなく「偽になる」例である。 -/
theorem bdDirection_breaks_prop_1_7 :
    ¬ ∀ (Pt : Type) (condE condD diffY : Pt → ℝ) (e : ℕ), 0 < e →
        (∀ x, condE x - condD x = diffY x) →
        (∀ x, 0 ≤ condE x) →
        (∀ x, 0 ≤ condD x) →
        BDle diffY (fun x => (1 - 1 / (e : ℝ)) * condE x) :=
  ABC3.Check.GenEll.prop_1_7_printed_direction_false

end ABC3.Gap.GenEll
