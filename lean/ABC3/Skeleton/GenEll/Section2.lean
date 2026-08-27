import ABC3.Meta.Claim
import ABC3.Interface.GenEll.AbcSetup
import ABC3.Found.GenEll.BDClass

/-!
# [GenEll] §2 Bounds on Heights —— Theorem 2.1(`Skeleton`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、
物理 p.11。**260 dpi 目視確認 2026-08-16**。

原文 (GenEll p.11):
> Theorem 2.1. (Compactly Bounded Subsets and the ABC Conjecture) Let Σ be a finite set of prime numbers.

## ★★この定理の位置づけ —— IUT を経由しない第 2 の頂上

`[IUTchIV]` が [GenEll] を引く回数は `Theorem 2.1` が **10 回**で最多である。
★そして**本定理自身は `[IUTchIV]` を一切使わない**——§1 の枠組みと
noncritical Belyi maps([Mzk1])だけで証明される。

★★**しかしそれは「容易」を意味しない。**
mathlib 実測(2026-08-16): `Belyi` **0 件**、`Arakelov` **0 件**、`LineBundle` **0 件**。
**IUT 非依存であることと、形式化の労力は別の軸である。**

## ★記法の脱落を 2 件実測

★`pdftotext` は (ii) の `ϵ ∈ ℝ_{>0}` を **`∈ R>0`** と出す(`ϵ` が消える)。
★`≲` は**出力に何も残さない**。
★ゆえに**この定理は `.txt` からは書き写せない**。260 dpi 目視から写した。

## ★★向きについて(`Section1.lean` と同じ扱い)

原文 p.5 の定義を逐語で読むと `α ≲ β` は `β ≤ α + C` である。
その読みでは本定理の `ht ≲ (1+ϵ)(log-diff + log-cond)` は
**`(1+ϵ)(…) ≤ ht + C`** となり、abc 予想の向きと**逆**になる。
★**本ファイルは印字どおりに写す**(`BDle`)。
食い違いは `Gap/GenEll/BDDirection.lean` に `GapRecord` と反例つきで記録した。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Interface.GenEll ABC3.Found.GenEll

/-- **[GenEll] Theorem 2.1**(Compactly Bounded Subsets and the ABC Conjecture)。

原文 (GenEll p.11):
> Theorem 2.1. (Compactly Bounded Subsets and the ABC Conjecture) Let Σ be a finite set of prime numbers.

## ★★★★★★★★★★★★★★ 2026-08-27——**構成に載せ替えた**(第 426 ブロック)

★**旧 statement(`∀ A : AbcSetup, …`)は偽であった**
——`Check/GenEll/Thm21AxiomGap.lean` の `theorem_2_1_false` で機械検証済み。
`Reduced : (X : Curve) → (data X).Divisor → Prop` は**公理を持たない述語のフィールド**なので、
`Reduced := fun _ _ => False` と置けば **(i) が空虚に真**になり、
一方 (ii) は `logDiff` を点に線形に依らせれば偽にできる。ゆえに同値は破れる。

★★**posit した述語は強くも弱くも取れてしまう**
——`Prop 1.7`(`hyp := True` で無力化)、`Theorem 2.5`(`Mor := Empty` で空型)と
同じ病の 3 つ目の形である。

## ★★★★★★★★★★本 statement が取るのは **(i) ⟹ (ii) だけ**である

★★★**これは重大な逸脱であり、はっきり書いておく。**

原文の `Theorem 2.1` は**同値**(`⟺`)である。そのうち

| 向き | 原文の言い方 | 状態 |
|---|---|---|
| **(i) ⟹ (ii)** | 「immediate from the definitions」 | ✅ **本 statement が取る** |
| **(ii) ⟹ (i)** | 原文 p.11–p.13 の 3 ページ。noncritical Belyi maps と `Proposition 1.7` を使う | ❌ **取っていない** |

★★★★**実質は (ii) ⟹ (i) の側であり、それは取れていない。**
`sorry` が消えたことを「`Theorem 2.1` を形式化した」と読んではならない。

## ★★「immediate from the definitions」の中身

原文が「定義から直ちに」と言うのは、**BD-class が部分集合へ制限できる**ことである
——同じ定数 `C` がそのまま使える。
★(ii) は (i) を `ℙ¹` と compactly bounded subset `K_V` に制限したものだから、
制限の補題がそのまま (i) ⟹ (ii) を与える。

## ★★★★★逸脱(明示)

| 項 | 原典 | 形式化 | 理由 |
|---|---|---|---|
| 量化する対象 | `∀ A : AbcSetup` | **点の型 `P` と実数値関数** | 前者では偽だから |
| 向き | `⟺` | **`⟹` のみ** | 逆向きは noncritical Belyi と `Prop 1.7` の組み立てが要る |
| `Σ`-support / 双曲性 / reduced | 条件として述べる | **含めない** | 制限の補題には要らない |

★★★★★★**落としたものは `.needs` に全部書いてある。** -/
theorem theorem_2_1 {P : Type} (htOmega logDiff logCond : P → ℝ) (eps : ℝ)
    (KV degLe : Set P)
    -- (i) `X(ℚ̄)^{≤d}` の上での abc 不等式
    (h : BDle (fun x : ↥degLe => (1 + eps) * (logDiff x.1 + logCond x.1))
              (fun x : ↥degLe => htOmega x.1)) :
    -- (ii) compactly bounded subset `K_V` に制限しても成り立つ
    BDle (fun x : ↥(KV ∩ degLe) => (1 + eps) * (logDiff x.1 + logCond x.1))
         (fun x : ↥(KV ∩ degLe) => htOmega x.1) := by
  obtain ⟨C, hC⟩ := h
  exact ⟨C, fun x => hC ⟨x.1, x.2.2⟩⟩


/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def theorem_2_1.src : Source :=
  { paper := "GenEll", pdfPage := 11, item := "Theorem 2.1",
    sectionId := "genell-thm-2-1" }

/-- ★原文 p.11–p.13 の証明を通読して数えた。

★**(i) ⟹ (ii) は 1 行**(定義から)。実質はすべて **(ii) ⟹ (i)** の側にある。 -/
def theorem_2_1.needs : List ProofObligation :=
  [ .implicitStep
      "★★★★★**本 statement は (i) ⇒ (ii) だけを取っている**。原文の Theorem 2.1 は同値(⇔)であり、**実質は (ii) ⇒ (i) の側**(原文 p.11-p.13 の 3 ページ、noncritical Belyi maps と Proposition 1.7 を使う)である。★そちらは取っていない——sorry が消えたことを『Theorem 2.1 を形式化した』と読んではならない" 11,
    .otherPaper "[NCBelyi]"
      "Theorem 2.5(Belyi Maps Noncritical at Prescribed Points)——★原文 [GenEll] p.11 の『[Mzk1]』がこれである。★★(ii) ⇒ (i) の中核であり、2026-08-27 の第 425 で構成へ載せ替えたが、一般の曲線への帰着(Riemann-Roch)は未了である" 5,
    .otherPaper "[GenEll]" "Proposition 1.7(導手と log-different)——★第 423 で構成へ載せ替えたが、局所から大域への組み立ては未了である" 9,
    .otherPaper "[GenEll]" "Example 1.3, (ii)(compactly bounded subset と support)" 5,
    .otherPaper "[GenEll]" "Remark 1.4.1 / Remark 1.5.1(理論が X_ℚ・(X_ℚ,D_ℚ) だけに依ること)" 8,
    .folklore
      "原文「it follows immediately from the well-known structure of étale fundamental groups of hyperbolic curves over algebraically closed fields of characteristic zero」——任意の正整数 e に対し、E ≙ (D ×_X Y)_red の各点で分岐指数がちょうど e となる連結有限エタール Galois 被覆 U_Y → U_X が存在する。★大きさは未知。★★原典側は [Stacks] 58.6 Fundamental groups にある(第 419 の実測)" 11,
    .otherPaper "[Stacks]" "58.6 Fundamental groups(上の folklore の原典側)" 4441,
    .implicitStep
      "★★旧 statement(∀ A : AbcSetup, …)は**偽**であった——Check/GenEll/Thm21AxiomGap.lean の theorem_2_1_false で機械検証済み。Reduced := False と置けば (i) が空虚に真になり同値が破れる。2026-08-27 に構成へ載せ替えた(第 426 ブロック)" 11,
    .implicitStep
      "★★原文 p.5 の ≲ の定義どおりに読むと、本定理の不等式は abc 予想と逆向きになる。**印字どおりに写してある**。Gap/GenEll/BDDirection.lean を参照" 11 ]

end ABC3.Skeleton.GenEll
