import ABC3.Interface.GenEll.GaloisRep
import Mathlib.Data.Real.Basic
import Mathlib.Data.Set.Finite.Basic

/-!
# [GenEll] §3–§4 大域理論 —— `M_ell(ℚ̄)` 上の高さの `Interface`

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、
物理 p.17–p.23。**260 dpi 目視確認 2026-08-16**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★`TorsionGaloisRepData` を**拡張する**(並立させない)

★`Interface/GenEll/GaloisRep.lean` の `TorsionGaloisRepData` は
`Theorem 3.8` の statement のためだけに作った**狭い**構造体である。
§3 の残りと §4 は同じ語(`M_ell(ℚ̄)`・`ht^Falt`・compactly bounded・Galois-finite)を
使うので、**別に posit すると 2 つの語彙が並立して腐る**。

★ゆえに `extends` で**機械的に接続する**。
`Theorem 3.8` の statement と `Corollary 4.3/4.4` の statement は、
**同じ `EllClass`・同じ `faltingsHeight` の上に立つことが型で保証される**。

## ★mathlib 実測(2026-08-16)

| 語 | mathlib |
|---|---|
| Faltings 高さ `ht^Falt` | ★**0 件** |
| `deg_∞` / `ht_∞`(無限遠因子) | ★**0 件** |
| `log-diff_{M̄_ell}` | ★**0 件**(`Definition 1.5, (iii)` を `M̄_ell` に適用したもの) |
| 楕円曲線のモジュライスタック `M̄_ell` | ★**0 件** |
-/

namespace ABC3.Interface.GenEll

open ABC3.Meta

/-- **[GenEll] §3–§4 の大域理論**を受ける `Interface`。

★`TorsionGaloisRepData`(`Theorem 3.8` 用)を**拡張**している。
追加したフィールドはすべて原文 p.17–p.23 に**実際に現れる語**である。 -/
structure EllModuliData extends TorsionGaloisRepData where
  /-- 原文 `deg_∞`(無限遠因子の次数、`Lemma 3.2, (ii)` の大域版)。 -/
  degInf : EllClass → ℝ
  /-- 原文 `ht_∞`(無限遠因子に付随する高さ)。 -/
  htInf : EllClass → ℝ
  /-- 原文 `log-diff_{M̄_ell}`(`Definition 1.5, (iii)` を `M̄_ell` に適用したもの)。 -/
  logDiffMell : EllClass → ℝ
  /-- 原文 `M_ell(ℚ̄)^{≤d}`(`Example 1.3, (i)`)。 -/
  degLe : ℕ → Set EllClass
  /-- 原文「`E_L` … with **semi-stable reduction** at all the finite primes of `L`」。 -/
  SemiStable : Curve → Prop
  /-- 原文「`E_L` admits an **l-cyclic** subgroup scheme `H_L ⊆ E_L`」。 -/
  HasLCyclic : Curve → ℕ → Prop
  /-- 原文「`L` is a **minimal field of definition** of the point `[E_L]`」。 -/
  MinimalField : Curve → Prop
  /-- 原文「The Galois representation `Gal(ℚ̄/L) → GL₂(ℤ_l)` associated to `E_L`
  is **surjective**」(`Corollary 4.3/4.4`, (b))。 -/
  ImageSurjective : Curve → ℕ → Prop
  /-- 原文「`l_•` is **prime** to the primes of `ℚ` that **ramify** in `L`, as well as to
  the **ramification indices** of primes of `ℚ` in `L`」(`Corollary 4.3/4.4`, (a))。 -/
  PrimeToRamification : Curve → ℕ → Prop
  /-- 原文「`E_L` has at least **one** prime of [bad] **multiplicative** reduction」。
  ★`TorsionGaloisRepData.HasPotMultRed`(**潜在的**乗法還元)とは**別の述語**である
  ——`Lemma 3.7` は両者を使い分けている(p.18 目視確認)。 -/
  HasMultRed : Curve → Prop
  /-- 原文「`l_∘`, `l_•` are **prime** to the **primes of potentially multiplicative
  reduction**」(`Corollary 4.3/4.4`, (a) の前半)。
  ★`PrimeToLocalHeights`(局所高さと素)とは**別の条件**である——
  原文 (a) は "as well as to" で 2 つを並べている(p.22 目視確認)。 -/
  PrimeToMultPrimes : Curve → ℕ → Prop

/-- ★Track B は何を作らねばならないか。 -/
def EllModuliData.waiting : WaitingFor :=
  { what := "楕円曲線のモジュライスタック M̄_ell と無限遠因子、Faltings 高さ ht^Falt、deg_∞ / ht_∞、log-diff_{M̄_ell}、および半安定還元・l-cyclic 部分群スキーム・最小定義体"
    trackB := "Found/GenEll — ★ht^Falt は Arakelov 理論(§1)を、SemiStable / HasLCyclic は Tate 曲線(Interface/GenEll/TateLocal.lean)を要求する。★log-diff は Definition 1.5, (iii) すなわち ADiv(F) と degNormalized の上に立ち、その 2 つは Found/GenEll/ArithDiv.lean に実装済みである" }

/-! ## ★出典の紐付け(`.src`) -/

def EllModuliData.src : Source :=
  { paper := "GenEll", pdfPage := 17, item := "Proposition 3.4",
    sectionId := "genell-prop-3-4" }

end ABC3.Interface.GenEll
