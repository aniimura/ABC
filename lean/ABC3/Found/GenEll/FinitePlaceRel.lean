import ABC3.Meta.Claim
import Mathlib.NumberTheory.NumberField.Basic
import Mathlib.RingTheory.DedekindDomain.AdicValuation
import Mathlib.RingTheory.Ideal.GoingUp
import Mathlib.NumberTheory.RamificationInertia.Basic

/-!
# [GenEll] §1 —— 有限素点の引き戻し(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> — where cv ∈ Z if v ∈ V(F )non and cv ∈ R if v ∈ V(F )arc [cf. [Szp], §1.1]. Here, if

## ★★`HeightOneSpectrum.comap` は mathlib に**無い**(2026-08-17 実測)

`IsDedekindDomain.HeightOneSpectrum` に**素点の引き戻し**は無い
(`DedekindDomain/` 配下を `comap` で検索し、`HeightOneSpectrum` の comap は 0 件)。

★無限素点の側には `NumberField.InfinitePlace.comap` が**ある**のに、
有限素点の側には無い——**非対称である**。

## ★何に要るか

`degNormalized` の**底変換不変性**の有限側。
`w : 𝕍(K)^non` の係数を `e(w|v)·n_v`(`v = w` の引き戻し)と定めるとき、
まず `w ↦ v` の写像が要る。

★アルキメデス側は `Found/GenEll/InfinitePlaceRel.lean` で済んでいる
(`sum_arc_base_change`)。

## ★構成は 3 行で済む

- `asIdeal` は `Ideal.comap`
- 素であることは `Ideal.comap_isPrime`
- ★**非零であることが唯一の中身**で、`Ideal.IsIntegral.comap_ne_bot`
  (整拡大では素イデアルの縮約が非零)が mathlib に**ある**
-/

namespace ABC3.Found.GenEll

open NumberField IsDedekindDomain

variable (F K : Type*) [Field F] [NumberField F] [Field K] [NumberField K] [Algebra F K]

/-- ★**有限素点の引き戻し** `𝕍(K)^non → 𝕍(F)^non`。

★mathlib に無いので作る(上の docstring)。
非零性は `Ideal.IsIntegral.comap_ne_bot`——**整拡大では素イデアルの縮約が非零**。 -/
def hosComap (w : HeightOneSpectrum (𝓞 K)) : HeightOneSpectrum (𝓞 F) where
  asIdeal := w.asIdeal.comap (algebraMap (𝓞 F) (𝓞 K))
  isPrime := Ideal.comap_isPrime _ _
  ne_bot := Ideal.IsIntegral.comap_ne_bot (𝓞 F) w.ne_bot

@[simp] theorem hosComap_asIdeal (w : HeightOneSpectrum (𝓞 K)) :
    (hosComap F K w).asIdeal = w.asIdeal.comap (algebraMap (𝓞 F) (𝓞 K)) := rfl

/-- ★`w` は自身の引き戻しの上にある(`v ≤ w ∩ 𝓞_F` の形)。 -/
theorem hosComap_le_comap (w : HeightOneSpectrum (𝓞 K)) :
    (hosComap F K w).asIdeal ≤ w.asIdeal.comap (algebraMap (𝓞 F) (𝓞 K)) := le_rfl

/-- ★引き戻した素点は `w` の下にある——`Ideal.LiesOver` の形。 -/
instance hosComap_liesOver (w : HeightOneSpectrum (𝓞 K)) :
    w.asIdeal.LiesOver (hosComap F K w).asIdeal :=
  ⟨rfl⟩

/-! ## ★出典の紐付け(`.src`) -/

def hosComap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4, item := "§1 地の文(算術因子 ADiv(F))",
    sectionId := "genell-adiv" }

end ABC3.Found.GenEll
