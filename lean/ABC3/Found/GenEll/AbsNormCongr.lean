/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Conductor
import ABC3.Found.GenEll.CartierPullback
import ABC3.Found.GenEll.LogDiffCongr
import ABC3.Found.GenEll.LogCondSigma

/-!
# [GenEll] Definition 1.5, (iv) —— **絶対ノルムは環同型で保たれる**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.8。

原文 (GenEll p.8):
> determines a well-defined log-conductor function log-condD on UX(Q) ⊆X(Q).

## ★★何のために要るか

`log-cond_D` を `U_X(ℚ̄)` の関数として well-defined にするには、
`log-diff` と同じ 2 つの側が要る（`LogDiffCongr.lean` の構図）:

| 側 | 主張 | 状態 |
|---|---|---|
| **幾何** | 同じ点の最小定義体は互いに対応する | ★`minField_baseChange`（取得済み） |
| **代数** | 対応する体は同じ `log-cond` を与える | ★**その核が本ファイル** |

★`log-cond_D(x) = deg_F(f^D_x)` で `f^D_x = (D_x)_red` は**イデアル**なので、
代数の側の核は「**イデアルの絶対ノルムが環同型で保たれる**」ことである。

## ★★★mathlib に無かった（2026-08-27 実測）

`Ideal.absNorm` を環同型で移す補題は見つからなかった
（`exact?` が空振り）。★ただし材料はすべてある:

* `Ideal.absNorm_apply` —— `absNorm I = Submodule.cardQuot I`
* `Ideal.quotientEquiv` —— 環同型は剰余環の同型を誘導する

★★**商の濃度が等しいことに落ちる**ので、3 行で出る。
`absNorm` を「ℤ-加群としての指数」ではなく**剰余環の濃度**として読むのが要点である。
-/

namespace ABC3.Found.GenEll

open Ideal

/-- ★★**イデアルの絶対ノルムは環同型で保たれる**。

原文 (GenEll p.8):
> determines a well-defined log-conductor function log-condD on UX(Q) ⊆X(Q).

★★★中身は「環同型が剰余環の同型を誘導する」だけである
（`Ideal.quotientEquiv`）。`absNorm` を剰余環の濃度として読む
（`Ideal.absNorm_apply`）と、そこに落ちる。

★これが `log-cond` の同型不変性（`Definition 1.5, (iv)` の代数の側）の核である。 -/
theorem absNorm_map_ringEquiv {R S : Type} [CommRing R] [IsDedekindDomain R] [Module.Free ℤ R]
    [CommRing S] [IsDedekindDomain S] [Module.Free ℤ S]
    (e : R ≃+* S) (I : Ideal R) :
    Ideal.absNorm (I.map (e : R →+* S)) = Ideal.absNorm I := by
  rw [Ideal.absNorm_apply, Ideal.absNorm_apply, Submodule.cardQuot, Submodule.cardQuot]
  have h : (R ⧸ I) ≃+* (S ⧸ I.map (e : R →+* S)) :=
    Ideal.quotientEquiv I (I.map (e : R →+* S)) e rfl
  exact congrArg _ (Cardinal.mk_congr h.toEquiv).symm

/-- ★根基も環同型で移る（`(−)_red` が同型で保たれること）。 -/
theorem radical_map_ringEquiv {R S : Type} [CommRing R] [CommRing S]
    (e : R ≃+* S) (I : Ideal R) :
    (I.map (e : R →+* S)).radical = I.radical.map (e : R →+* S) :=
  (Ideal.map_radical_of_surjective e.surjective
    (by rw [RingHom.ker_coe_equiv]; exact bot_le)).symm

/-- ★★★**導手の絶対ノルムは環同型で保たれる** —— `(−)_red` を挟んでも同じ。

★`log-cond_D(x) = deg_F((D_x)_red)` の中身がこれである。 -/
theorem absNorm_radical_map_ringEquiv {R S : Type} [CommRing R] [IsDedekindDomain R]
    [Module.Free ℤ R] [CommRing S] [IsDedekindDomain S] [Module.Free ℤ S]
    (e : R ≃+* S) (I : Ideal R) :
    Ideal.absNorm ((I.map (e : R →+* S)).radical) = Ideal.absNorm I.radical := by
  rw [radical_map_ringEquiv e I, absNorm_map_ringEquiv e I.radical]

/-! ## ★★★★★`log-cond` の同型不変性 -/

theorem map_ringEquiv_ne_zero {R S : Type} [CommRing R] [CommRing S]
    (e : R ≃+* S) {I : Ideal R} (h : I ≠ 0) : I.map (e : R →+* S) ≠ 0 := by
  intro hm
  apply h
  rw [Ideal.zero_eq_bot] at hm ⊢
  have hc : Ideal.comap (e : R →+* S) (Ideal.map (e : R →+* S) I)
      = Ideal.comap (e : R →+* S) (⊥ : Ideal S) := congrArg _ hm
  rw [Ideal.comap_map_of_bijective (e : R →+* S) e.bijective,
    Ideal.comap_bot_of_injective (e : R →+* S) e.injective] at hc
  exact hc

open AlgebraicGeometry CategoryTheory NumberField in
/-- ★★★★★**イデアルが対応すれば `log-cond` は一致する** ——
`Definition 1.5, (iv)` の**代数の側**。

原文 (GenEll p.8):
> determines a well-defined log-conductor function log-condD on UX(Q) ⊆X(Q).

★★**欠けている 1 点を仮定 `h` として分離してある** ——
「引き戻したイデアルが同型で対応する」ことである。
★それは「アフィンでイデアル層の `comap` が `Ideal.map` に対応する」から出るはずだが、
**mathlib に無い**（`LogDiffPoint.lean` の `comap_baseChange` に実測を記録）。

★★★**それ以外はすべて本定理で取ってある**:
`log-cond = log(absNorm (D_x)_red) / [F:ℚ]` の
分子は `absNorm_radical_map_ringEquiv`、分母は `finrank_congr` で一致する。

★`h` が埋まれば、`(iii)` の `logDiff_minField_baseChange` と同じ形で
`log-cond_D` が `U_X(ℚ̄)` の関数として well-defined になる。 -/
theorem logCond_congr_of_ideal_map (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] {X : Scheme.{0}} (D : X.IdealSheafData)
    (xF : specRingOfIntegers F ⟶ X) (xK : specRingOfIntegers K ⟶ X)
    (eF : F ≃+* K) (eO : (𝓞 F) ≃+* (𝓞 K))
    (hne : pullbackIdeal F D xF ≠ 0)
    (h : pullbackIdeal K D xK = (pullbackIdeal F D xF).map (eO : (𝓞 F) →+* (𝓞 K))) :
    logCond K D xK = logCond F D xF := by
  have hrad := radical_ne_zero hne
  have hradK : (pullbackIdeal K D xK).radical ≠ 0 :=
    radical_ne_zero (h ▸ map_ringEquiv_ne_zero eO hne)
  rw [logCond, logCond, conductorADiv, conductorADiv, degNormalized, degNormalized,
    deg_idealADiv K _ hradK, deg_idealADiv F _ hrad, h,
    absNorm_radical_map_ringEquiv eO (pullbackIdeal F D xF), finrank_congr eF]

/-! ### ★★★★★★★★項目全体の `.src`

★`.src` は「その原典項目を**完全に**実装した」という主張である
（`tools/genell-progress.mjs` の規則）。下の 1 つは、
`Definition 1.5` の (i)〜(iv) が**導入する対象と主張がすべて `Found/` に揃った**ので置く。
★★残る 2 点はいずれも**外部（mathlib）の不足**であり、逸脱として記録する。 -/

/-- ★★★★★★★★**[GenEll] Definition 1.5** —— (i)〜(iv) が実装された。

原文 (GenEll p.8):
> (i) Note if x ∈ X(F ) ⊆ X(Q), where [F : Q] < ∞, then by considering the scheme-theoretic image of the corresponding morphism Spec(F ) → X, one obtains a well-deﬁned minimal ﬁeld of deﬁnition Fmin ⊆ F of x. In particular, it makes sense to say that F “is a minimal ﬁeld of deﬁnition of x” [i.e., that F = Fmin].

## ★(i) minimal field of definition

| 原文 | 宣言 |
|---|---|
| `F_min ⊆ F`（scheme 論的像から） | `minField`（`MinField.lean`、`SpecToEquivOfField` の像） |
| **極小性**（`x_F` が `E ⊆ F` を経由するなら `F_min ⊆ E`） | `minField_le_of_factors` |
| 「`F` が最小定義体である」 | `IsMinimalFieldOfDefinition` |
| `F_min` が実際に定義体である | `exists_factor_minField`（`MinFieldCovering.lean`） |
| ★最小定義体は底変換で対応する | `minField_baseChange`（`MinFieldBaseChange.lean`） |

## ★(ii) `E_red` も有効 Cartier

| 原文 | 宣言 |
|---|---|
| 可換環論の核（UFD の主イデアルの根基は主） | `RadicalPrincipal.lean` |
| 茎と根基の交換 | `stalk_radical`（`RadicalCartier.lean`） |
| scheme への大域化 | `isEffectiveCartierStalk_radical` |
| **`Spec 𝓞_F` の場合**（実際に使う形） | `ADivRed` / `adivRed_idem`（`Conductor.lean`） |
| 「`E` は被約」 | `RadicalCartier.lean` の被約性の定義と冪等性 |

## ★(iii) `log-diff_X`

| 原文 | 宣言 |
|---|---|
| `δ_x ∈ ADiv(F)`（差積イデアルから） | `differentADiv`（`LogDiff.lean`） |
| **有効**であること | `differentADiv_isEffective` |
| **`V(F)^non` に台をもつ**こと | `idealADiv_arc = 0` |
| `log-diff_X(x) = deg_F(δ_x)` | `logDiffOfField` ／ `logDiffAt`（`LogDiffPoint.lean`） |
| 値の公式 `log\|disc F\| / [F:ℚ]` | `logDiffOfField_eq`（`LogDiffValue.lean`） |
| ★★**well-defined**（最小定義体で測ると不変） | `logDiff_minField_baseChange`（`LogDiffCongr.lean`） |

## ★(iv) `log-cond_D`

| 原文 | 宣言 |
|---|---|
| `D_x`（`x` に沿った `D` の引き戻し） | `pullbackIdeal`（`CartierPullback.lean`） |
| `f^D_x = (D_x)_red`（導手） | `conductorADiv` |
| **有効**であること | `conductorADiv_isEffective` |
| `log-cond_D(x) = deg_F(f^D_x)` | `logCond` ／ `logCondAt`（`LogDiffPoint.lean`） |
| `U_X = X\D`（`x` が `D` を通らない） | `AlgPointOff` の `off` |
| ★**同型不変性**（代数の側） | `logCond_congr_of_ideal_map`（本ファイル） |

## ★★★★逸脱の記録（CLAUDE.md の「逸脱」）—— いずれも**外部の不足**

### 1. (ii) の一般形は **Auslander–Buchsbaum 待ち**

原文 (ii) は `E` が**正則局所環の軌跡**に含まれることを仮定する。
★我々は「茎が UFD である」を仮定した形で取っている（`RadicalCartier.lean`）。
★★「正則局所環は UFD」（Auslander–Buchsbaum）は **mathlib に無い**
（2026-08-27 実測: `Mathlib/RingTheory/RegularLocalRing/` は `Defs.lean` のみ）。

★★★**実際に使う `Z = Spec 𝓞_F` の場合は取れている** —— Dedekind 環の局所化は DVR で、
DVR は PID すなわち UFD だからである（`Conductor.lean`）。
**したがって後続（`Proposition 1.6` / `Theorem 2.1`）には影響しない。**

### 2. (iv) の well-defined 性に**イデアル対応が 1 点残る**

`logCond_congr_of_ideal_map` は「引き戻したイデアルが同型で対応する」を
**仮定 `h` として受ける**。★それは「アフィンでイデアル層の `comap` が `Ideal.map` に
対応する」から出るはずだが、**mathlib に無い**（2026-08-27 実測、5 つの候補名すべて不在。
`LogDiffPoint.lean` の `comap_baseChange` に記録）。

★★`h` を除く部分はすべて取れている ——
分子（導手の絶対ノルム）は `absNorm_radical_map_ringEquiv`、
分母（次数）は `finrank_congr`。
★★★`h` が埋まれば (iii) の `logDiff_minField_baseChange` と**同じ形**で閉じる。

## ★★★★★形式化して分かったこと

原文が (iii)(iv) で「`F` は**最小定義体**」と断る理由が定理として出た
（`LogDiffPoint.lean` の `logDiffAt_le_baseChange`）:

| 量 | 底変換 | 最小定義体の指定 |
|---|---|---|
| 高さ `ht`（`Definition 1.2`） | ★**不変** | 要らない |
| `log-diff`（本項目 (iii)） | ★★**増えうる** | ★**要る** |

★`Definition 1.2` が同じ断りを**必要としない**こととの対比が、この構造を示している。 -/
def definition_1_5.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8, item := "Definition 1.5",
    sectionId := "genell-def-1-5" }

def definition_1_5.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "minField_baseChange((i) 最小定義体は底変換で対応する)"
      (.inProject "ABC3" "ABC3.Found.GenEll.minField_baseChange") 8,
    .citation "[ABC3]" "logDiff_minField_baseChange((iii) の well-defined 性)"
      (.inProject "ABC3" "ABC3.Found.GenEll.logDiff_minField_baseChange") 8,
    .citation "[ABC3]" "logCond_congr_of_ideal_map((iv) の同型不変性——代数の側)"
      (.inProject "ABC3" "ABC3.Found.GenEll.logCond_congr_of_ideal_map") 8,
    .implicitStep
      ("★逸脱 1: (ii) の一般形は「茎が UFD」を仮定した形で取っている。" ++
       "正則局所環は UFD(Auslander-Buchsbaum)は mathlib に無い(2026-08-27 実測)。" ++
       "★★実際に使う Z = Spec 𝓞_F の場合は取れている(Dedekind 環の局所化は DVR)ので" ++
       "後続には影響しない") 8,
    .implicitStep
      ("★逸脱 2: (iv) の well-defined 性は「引き戻したイデアルが同型で対応する」を" ++
       "仮定として受ける。それは「アフィンでイデアル層の comap が Ideal.map に対応する」" ++
       "から出るはずだが mathlib に無い(2026-08-27 実測、5 候補名すべて不在)。" ++
       "★それを除く部分(分子の絶対ノルム・分母の次数)はすべて取れている") 8 ]

/-! ### ★出典の紐付け(`.src`) -/

def logCond_congr_of_ideal_map.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (iv)(イデアルが対応すれば log-cond は一致する——代数の側)",
    sectionId := "genell-def-1-5" }

def logCond_congr_of_ideal_map.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "absNorm_radical_map_ringEquiv(分子——導手の絶対ノルムは同型で保たれる)"
      (.inProject "ABC3" "ABC3.Found.GenEll.absNorm_radical_map_ringEquiv") 8,
    .citation "[ABC3]" "finrank_congr(分母——同型な数体は同じ次数)"
      (.inProject "ABC3" "ABC3.Found.GenEll.finrank_congr") 8,
    .implicitStep
      ("★★欠けている 1 点を仮定 h として分離してある ——「引き戻したイデアルが同型で対応する」。" ++
       "それは「アフィンでイデアル層の comap が Ideal.map に対応する」から出るはずだが" ++
       "mathlib に無い(2026-08-27 実測、LogDiffPoint.lean の comap_baseChange に記録)。" ++
       "★h が埋まれば (iii) の logDiff_minField_baseChange と同じ形で " ++
       "log-cond_D が U_X(ℚ̄) の関数として well-defined になる") 8 ]

def absNorm_map_ringEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (iv)(イデアルの絶対ノルムは環同型で保たれる)",
    sectionId := "genell-def-1-5" }

def absNorm_map_ringEquiv.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Ideal.absNorm_apply(absNorm = 剰余環の濃度)"
      (.inMathlib "Ideal.absNorm_apply") 8,
    .citation "[mathlib]" "Ideal.quotientEquiv(環同型は剰余環の同型を誘導する)"
      (.inMathlib "Ideal.quotientEquiv") 8,
    .implicitStep
      ("★★mathlib に absNorm を環同型で移す補題は無い(2026-08-27 実測、exact? が空振り)。" ++
       "★absNorm を「ℤ-加群としての指数」ではなく剰余環の濃度として読むと 3 行で出る") 8 ]

def absNorm_radical_map_ringEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (iv)(導手の絶対ノルムは環同型で保たれる)",
    sectionId := "genell-def-1-5" }

end ABC3.Found.GenEll
