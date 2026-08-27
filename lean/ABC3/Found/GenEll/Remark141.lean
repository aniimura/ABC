/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.PointCorrespondence
import ABC3.Found.GenEll.HeightClass
import ABC3.Found.GenEll.VerticalTwist
import ABC3.Found.GenEll.VerticalBound
import ABC3.Found.GenEll.IdealComparable

/-!
# [GenEll] Remark 1.4.1 —— **2 つの `ℤ`-モデルにまたがる形**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.8。

原文 (GenEll p.8):
> Remark 1.4.1. Observe that it follows immediately from the definitions, together with Proposition 1.4, (iii), that the theory of

## ★★★★★★★★何を足したか

`HeightClass.lean` の `htArith_bdeq_of_pullbackAPic_eq` は
**同じ `X` の上の 2 つの `ArithCartier`** について述べていた。
★`Skeleton/GenEll/Section1.lean` の `remark_1_4_1` もその形で、
逸脱として「異なる 2 つの `ℤ`-モデルの比較は含めていない」と明記していた。

★★`PointCorrespondence.lean` で**点の対応 `ePt` が構成できる**ようになったので、
**2 つの `ℤ`-モデルにまたがる形**が書ける:

    ∀ x,  x_F^*D の類 = (ePt x)_F^*E の類   ⟹   ht_D ≈ ht_E ∘ ePt

★★★定数は `0`——類が一致すれば高さは**そのもの**が一致する。

## ★★★★残っている段（明示）

★原文の「`X_ℚ` だけに依る」を**仮定なしで**出すには、
「2 つのモデルの算術因子が `ℚ` 上対応すれば、引き戻した類も対応する」が要る。
★★これは `ℤ[1/Σ]` の上では `IsoDescent`／`DivisorDescent` で出るが、
`Σ` の上（＝**垂直因子**の差）が残る。

★★★**その差は一様に有界である**——垂直因子は `Spec ℤ` 上の因子の引き戻しなので、
`x^*(X_q) = (q)` となり `deg_F` は `log q`（点にも定義体にも依らない）。

★★★★★**この段は `VerticalTwist.lean` で取れた**（本セッション、第 359 ブロック）:
差が `Spec ℤ` 上の `N : ADiv ℚ` の底変換なら、高さの差は**定数 `deg_ℚ(N)`**。
原文が `Proposition 1.4, (ii)` の証明で 1 行で使っている段そのものである。

★残るのは「**2 つのモデルのずれが必ずその形になる**」——
ファイバーが既約なら垂直因子は `∑_q n_q · f^*(q)` の形になるが、
可約な場合は交点数の評価が要る。★本ファイルはそこを仮説に置いている。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Limits

variable (F : Type) [Field F] [NumberField F]

/-- ★★★★★★★★**高さは 2 つの `ℤ`-モデルにまたがっても、引き戻した類だけで決まる**。

原文 (GenEll p.8):
> Remark 1.4.1. Observe that it follows immediately from the definitions, together with Proposition 1.4, (iii), that the theory of

★`HeightClass.lean` の同名補題（同じ `X` の上）を、**点の対応 `ePt` を挟んで**
2 つのモデルへ広げたもの。

★★定数は `0`——類が一致すれば高さはそのものが一致する。
★★★これが `Proposition 1.4, (iii)`（高さの BD-class は類だけに依る）の
`Remark 1.4.1` 版である。 -/
theorem htArith_bdeq_of_pullbackAPic_eq' {X X' : Scheme.{0}}
    (D : ArithCartier X) (E : ArithCartier X')
    (ePt : (specRingOfIntegers F ⟶ X) → (specRingOfIntegers F ⟶ X'))
    (h : ∀ xF, pullbackAPic F D xF = pullbackAPic F E (ePt xF)) :
    BDeq (fun xF => htArith F D xF) (fun xF => htArith F E (ePt xF)) :=
  ⟨0, fun xF => by
    show |htArith F D xF - htArith F E (ePt xF)| ≤ 0
    rw [htArith_eq_degNormalizedAPic, htArith_eq_degNormalizedAPic, h xF, sub_self, abs_zero]⟩

/-- ★★★★★★★★★**[GenEll] Remark 1.4.1 —— 点の対応まで込みの形**。

原文 (GenEll p.8):
> Remark 1.4.1. Observe that it follows immediately from the definitions, together with Proposition 1.4, (iii), that the theory of

段 `n` での同型 `φ` から**点の対応 `ePt` を構成し**、
その上で「引き戻した類が一致すれば高さは BD-同値」を述べる。

★★`ePt` はもう仮定ではない（`PointCorrespondence.lean` の `exists_ePt`）。
★★★残る仮定 `h`（類が一致する）を落とすには**垂直因子の差の評価**が要る
——ファイル冒頭の「残っている段」を見よ。 -/
theorem remark_1_4_1_of_descent {X X' : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (f' : X' ⟶ Spec (CommRingCat.of ℤ)) [IsProper f']
    {n : ℕᵒᵖ}
    (A : Type) [CommRing A] [IsDomain A] [Algebra (NumberField.RingOfIntegers F) A]
    [Algebra A F] [IsScalarTower (NumberField.RingOfIntegers F) A F] [IsFractionRing A F]
    (hinv : IsUnit (algebraMap ℤ A ((Nat.factorial n.unop : ℕ) : ℤ)))
    (φ : bcObj f n ⟶ bcObj f' n)
    (D : ArithCartier X) (E : ArithCartier X') :
    ∃ ePt : (specRingOfIntegers F ⟶ X) → (specRingOfIntegers F ⟶ X'),
      (∀ xF, Spec.map (CommRingCat.ofHom (algebraMap (NumberField.RingOfIntegers F) A)) ≫ ePt xF
        = liftPointToBc F A hinv f xF ≫ φ ≫
          pullback.snd (overRatTowerDiagram.obj n).hom f')
      ∧ ((∀ xF, pullbackAPic F D xF = pullbackAPic F E (ePt xF)) →
        BDeq (fun xF => htArith F D xF) (fun xF => htArith F E (ePt xF))) := by
  obtain ⟨ePt, hcompat⟩ := exists_ePt F f f' A hinv φ
  exact ⟨ePt, hcompat, fun h => htArith_bdeq_of_pullbackAPic_eq' F D E ePt h⟩

/-- ★★★★★★★★★★**[GenEll] Remark 1.4.1 —— 垂直なひねりを許した形**。

原文 (GenEll p.8):
> Remark 1.4.1. Observe that it follows immediately from the definitions, together with Proposition 1.4, (iii), that the theory of

上の `remark_1_4_1_of_descent` は「引き戻した**類が一致する**」を仮定していた。
★★**その仮定を弱められる**——`VerticalTwist.lean` により、
差が `Spec ℤ` から来る固定の `N : ADiv ℚ` の底変換であれば、
高さの差は**定数 `deg_ℚ(N)`** になる。

★★★これは原文が `Proposition 1.4, (ii)` の証明で 1 行で使っている段
（`ht_{L⊗M} ≈ ht_L`、`M` は `Spec ℤ` から来る）そのものである。

★★★★**2 つの `ℤ`-モデルのずれは、まさにこの「垂直な」方向にしか起きない**
——生成ファイバーが同型だからである。 -/
theorem remark_1_4_1_of_descent_twist {X X' : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (f' : X' ⟶ Spec (CommRingCat.of ℤ)) [IsProper f']
    {n : ℕᵒᵖ}
    (A : Type) [CommRing A] [IsDomain A] [Algebra (NumberField.RingOfIntegers F) A]
    [Algebra A F] [IsScalarTower (NumberField.RingOfIntegers F) A F] [IsFractionRing A F]
    (hinv : IsUnit (algebraMap ℤ A ((Nat.factorial n.unop : ℕ) : ℤ)))
    (φ : bcObj f n ⟶ bcObj f' n)
    (D : ArithCartier X) (E : ArithCartier X') (N : ADiv ℚ) :
    ∃ ePt : (specRingOfIntegers F ⟶ X) → (specRingOfIntegers F ⟶ X'),
      (∀ xF, Spec.map (CommRingCat.ofHom (algebraMap (NumberField.RingOfIntegers F) A)) ≫ ePt xF
        = liftPointToBc F A hinv f xF ≫ φ ≫
          pullback.snd (overRatTowerDiagram.obj n).hom f')
      ∧ ((∀ xF, pullbackADiv F D xF - pullbackADiv F E (ePt xF) = baseChange ℚ F N) →
        BDeq (fun xF => htArith F D xF) (fun xF => htArith F E (ePt xF))) := by
  obtain ⟨ePt, hcompat⟩ := exists_ePt F f f' A hinv φ
  exact ⟨ePt, hcompat, fun h => htArith_bdeq_of_baseChange' F D E ePt N h⟩

/-- ★★★★★★★★★★**[GenEll] Remark 1.4.1 —— 到達点**。

原文 (GenEll p.8):
> Remark 1.4.1. Observe that it follows immediately from the definitions, together with Proposition 1.4, (iii), that the theory of

仮定は 3 段で弱まった:

| 版 | 仮定 | 定数 |
|---|---|---|
| `remark_1_4_1_of_descent` | 引き戻した**類が一致**する | `0` |
| `remark_1_4_1_of_descent_twist` | 差が `Spec ℤ` から来る（`VerticalTwist`） | `deg_ℚ(N)` |
| **本定理** | 差のイデアルが**有理整数 `n` を含む**（`VerticalBound`） | `log n` |

★★★★★最後の版が**幾何が実際に与える形**である——
垂直因子 `V ≤ m·X_q` の引き戻しは `q^m` を含むので、
**ファイバーが可約でも交点数は要らない**。

★★定数 `log n` は**点にも定義体にも依らない**。 -/
theorem remark_1_4_1_of_descent_bound {X X' : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (f' : X' ⟶ Spec (CommRingCat.of ℤ)) [IsProper f']
    {m : ℕᵒᵖ}
    (A : Type) [CommRing A] [IsDomain A] [Algebra (NumberField.RingOfIntegers F) A]
    [Algebra A F] [IsScalarTower (NumberField.RingOfIntegers F) A F] [IsFractionRing A F]
    (hinv : IsUnit (algebraMap ℤ A ((Nat.factorial m.unop : ℕ) : ℤ)))
    (φ : bcObj f m ⟶ bcObj f' m)
    (D : ArithCartier X) (E : ArithCartier X') (n : ℕ) (hn : n ≠ 0) :
    ∃ ePt : (specRingOfIntegers F ⟶ X) → (specRingOfIntegers F ⟶ X'),
      (∀ xF, Spec.map (CommRingCat.ofHom (algebraMap (NumberField.RingOfIntegers F) A)) ≫ ePt xF
        = liftPointToBc F A hinv f xF ≫ φ ≫
          pullback.snd (overRatTowerDiagram.obj m).hom f')
      ∧ (∀ J : (specRingOfIntegers F ⟶ X) → Ideal (NumberField.RingOfIntegers F),
          (∀ xF, ((n : ℕ) : NumberField.RingOfIntegers F) ∈ J xF) →
          (∀ xF, pullbackADiv F E (ePt xF) - pullbackADiv F D xF = idealADiv F (J xF)) →
          BDeq (fun xF => htArith F D xF) (fun xF => htArith F E (ePt xF))) := by
  obtain ⟨ePt, hcompat⟩ := exists_ePt F f f' A hinv φ
  exact ⟨ePt, hcompat, fun J hJn h => htArith_bdeq_of_idealADiv_diff F D E ePt n hn J hJn h⟩

/-- ★★★★★★★★★★**[GenEll] Remark 1.4.1 —— 最も弱い仮定の版**。

原文 (GenEll p.8):
> Remark 1.4.1. Observe that it follows immediately from the definitions, together with Proposition 1.4, (iii), that the theory of

仮定の梯子（弱いほど下）:

| 版 | 仮定 | 定数 | 出どころ |
|---|---|---|---|
| `remark_1_4_1_of_descent` | 引き戻した**類が一致** | `0` | `HeightClass` |
| `remark_1_4_1_of_descent_twist` | 差が `Spec ℤ` から来る | `deg_ℚ(N)` | `VerticalTwist` |
| `remark_1_4_1_of_descent_bound` | 差のイデアルが `n` を含む | `log n` | `VerticalBound` |
| **本定理** | 2 つのイデアルが **`n` で互いに比較できる** | `log n` | `IdealComparable` |

★★★★★最下段が**「`n` を反転すれば一致する」の定量版**である
——これが降下（`exists_pair_descent`）が実際に与える形に最も近い。

★★★★★★**一様な `n` をスキームから取る段は閉じた**（2026-08-27、第 365 ブロック）:
`ChartwiseUniform.lean` の `htArith_bdeq_of_chartwise_localization` が

> `X` が有限個のネーターなアフィンチャートで覆われ、
> 各チャートで `D` と `E` が `N` を反転して一致するなら `ht_D ≈ ht_E`

を与える。★`mathlib-gap.json` の `vertical-divisors-generate-pic-kernel` は
**5/5 で閉じた**。

★★残るのは原文**第 2 文**の側だけである
——「this theory may be applied to any normal, projective scheme `Y` over `ℚ`」は
射影埋め込みから `ℤ`-モデルを作る段を要求し、それは
`ample-and-projective-embedding` の在庫不足である。 -/
theorem remark_1_4_1_of_descent_comparable {X X' : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (f' : X' ⟶ Spec (CommRingCat.of ℤ)) [IsProper f']
    {m : ℕᵒᵖ}
    (A : Type) [CommRing A] [IsDomain A] [Algebra (NumberField.RingOfIntegers F) A]
    [Algebra A F] [IsScalarTower (NumberField.RingOfIntegers F) A F] [IsFractionRing A F]
    (hinv : IsUnit (algebraMap ℤ A ((Nat.factorial m.unop : ℕ) : ℤ)))
    (φ : bcObj f m ⟶ bcObj f' m)
    (D : ArithCartier X) (E : ArithCartier X') (n : ℕ) (hn : n ≠ 0) :
    ∃ ePt : (specRingOfIntegers F ⟶ X) → (specRingOfIntegers F ⟶ X'),
      (∀ xF, Spec.map (CommRingCat.ofHom (algebraMap (NumberField.RingOfIntegers F) A)) ≫ ePt xF
        = liftPointToBc F A hinv f xF ≫ φ ≫
          pullback.snd (overRatTowerDiagram.obj m).hom f')
      ∧ ((∀ xF, pullbackIdeal F D.divisor xF ≠ 0) →
         (∀ xF, pullbackIdeal F E.divisor (ePt xF) ≠ 0) →
         (∀ xF, Ideal.span {((n : ℕ) : NumberField.RingOfIntegers F)}
             * pullbackIdeal F D.divisor xF ≤ pullbackIdeal F E.divisor (ePt xF)) →
         (∀ xF, Ideal.span {((n : ℕ) : NumberField.RingOfIntegers F)}
             * pullbackIdeal F E.divisor (ePt xF) ≤ pullbackIdeal F D.divisor xF) →
         (∀ xF, (archADiv F D.green xF).sum (fun _ r => r)
             = (archADiv F E.green (ePt xF)).sum (fun _ r => r)) →
         BDeq (fun xF => htArith F D xF) (fun xF => htArith F E (ePt xF))) := by
  obtain ⟨ePt, hcompat⟩ := exists_ePt F f f' A hinv φ
  exact ⟨ePt, hcompat, fun hD0 hE0 h1 h2 harc =>
    htArith_bdeq_of_ideal_comparable F D E ePt n hn hD0 hE0 h1 h2 harc⟩

/-! ### ★出典の紐付け(`.src`)

★★**条つきである。** 原文の「`X_ℚ` だけに依る」を仮定なしで出すには
一様な `n` をスキームの側から取る段が要る。 -/

def htArith_bdeq_of_pullbackAPic_eq'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Remark 1.4.1(高さは 2 つの ℤ-モデルにまたがっても引き戻した類だけで決まる)",
    sectionId := "genell-rem-1-4-1" }

def remark_1_4_1_of_descent.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Remark 1.4.1(点の対応を構成した形。垂直因子の差の評価は含まない)",
    sectionId := "genell-rem-1-4-1" }

def remark_1_4_1_of_descent_twist.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Remark 1.4.1(垂直なひねりを許した形。差が Spec ℤ から来る場合)",
    sectionId := "genell-rem-1-4-1" }

def remark_1_4_1_of_descent_bound.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Remark 1.4.1(到達点——差のイデアルが有理整数 n を含む場合、定数 log n)",
    sectionId := "genell-rem-1-4-1" }

def remark_1_4_1_of_descent_comparable.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Remark 1.4.1(最弱形——引き戻したイデアルが n で互いに比較できる場合)",
    sectionId := "genell-rem-1-4-1" }

def remark_1_4_1_of_descent.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_ePt(点の対応を固有性から構成する)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_ePt") 8,
    .citation "[ABC3]" "htArith_eq_degNormalizedAPic(高さは類の正規化次数)"
      (.inProject "ABC3" "ABC3.Found.GenEll.htArith_eq_degNormalizedAPic") 6,
    .otherPaper "[GenEll]" "Proposition 1.4, (iii)(高さの BD-class は L_ℚ の同型類だけに依る)" 6,
    .implicitStep
      ("★★★残る段: 原文の「X_ℚ だけに依る」を仮定なしで出すには" ++
       "「2 つのモデルの算術因子が ℚ 上対応すれば引き戻した類も対応する」が要る。" ++
       "★ℤ[1/Σ] の上では IsoDescent／DivisorDescent で出るが、Σ の上" ++
       "(＝垂直因子の差)が残る") 8,
    .citation "[ABC3]" "htArith_bdeq_of_baseChange′(垂直なひねりは高さを定数だけ動かす)"
      (.inProject "ABC3" "ABC3.Found.GenEll.htArith_bdeq_of_baseChange'") 6,
    .citation "[ABC3]"
      "htArith_bdeq_of_chartwise_localization(垂直な差の一様有界性——完全な形)"
      (.inProject "ABC3" "ABC3.Found.GenEll.htArith_bdeq_of_chartwise_localization") 6,
    .implicitStep
      ("★★★★★★**垂直な差の一様有界性は 2026-08-27 に閉じた**" ++
       "(mathlib-gap の vertical-divisors-generate-pic-kernel、5/5)。" ++
       "★当初「Cl(X) → Cl(U) の完全列が要る」と登記したのは過大であった" ++
       "——実際に使ったのは 有限生成性・準コンパクト性・absNorm の単調性・" ++
       "Spec ℤ が終対象であること の 4 つだけである。" ++
       "★★**交点数も Weil 因子類群も要らなかった**") 8,
    .implicitStep
      ("★★原文の第 2 文「In particular, this theory may be applied to any normal, " ++
       "projective scheme Y over ℚ」は、射影埋め込みから ℤ-モデルを作る段" ++
       "(P^n_ℤ での閉包)を要求する。★mathlib-gap の ample-and-projective-embedding" ++
       "と同じ在庫不足である") 8 ]

end ABC3.Found.GenEll
