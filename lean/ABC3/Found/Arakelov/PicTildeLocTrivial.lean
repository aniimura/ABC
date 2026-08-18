import ABC3.Found.Arakelov.PicResGen
import ABC3.Found.Arakelov.PicPointwise

/-!
# Arakelov (B1) 第 132 ブロック —— ★★★★★★★★**`tilde M` は局所自明である**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★可逆加群の層は局所自明である

第 93 ブロックから 40 ブロックかけて積み上げた道具立てが、ここで**閉じる**。

    可逆加群 M  ⟹  tilde M は 𝒪_{Spec R} 加群として局所自明

## ★★組み上げの筋

| 段 | ブロック | 内容 |
|---|---|---|
| 1 | 128 | 点ごとに調べれば良い(`isLocallyTrivial_of_pointwise`) |
| 2 | ——  | 各点 `x` の近傍に `M_g ≅ R_g` となる `D(g)`(`exists_away_linearEquiv`) |
| 3 | 127 | `D(g)` の上の切断 `s₀` が生成元 |
| 4 | 122 | `D(h·g)` たちは `Over (D g)` の site を覆う |
| 5 | 120 | 制限と局所化の可換図式 |
| 6 | 130 | ★★**制限した生成元は生成元**(`M` の言葉) |
| 7 | 131 | ★★同上(切断の言葉) |
| 8 | 115 | 生成族の上で全単射なら同型 |

## ★★★型の 2 経路との闘い(記録)

★組み上げで詰まった点は**すべて型の経路**であった([[ring-instance-two-paths]]):

| 症状 | 逃げ道 |
|---|---|
| `rw [hgen]` が motive で落ちる | ★`hgen` の**左辺を `s0` にする**(型が合う側から書く) |
| `Opens (PrimeSpectrum R)` vs `(Spec R).Opens` | ★`(X := Spec R)` を明示 |
| `Unique` の instance が見つからない | ★★`@restrictSec … (overTermInst _ _)` で**手渡し** |

★★★どれも「`exact`/`@` で項を手渡せば通る」——[[exact-term-over-rw]] の 4 例目である。

## ★★★★これで `PicardData` の最後の欄が埋まる

`IsInvertibleSheaf (tilde M)` の第 2 項(局所自明性)がこれである。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

noncomputable section

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★**可逆加群の層は局所自明である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★40 ブロック(93–131)の道具立てがここで閉じる。 -/
theorem isLocallyTrivial_tilde (R : CommRingCat.{u}) (M : ModuleCat.{u} (R : Type u))
    [Module.Invertible (R : Type u) (M : Type u)] :
    IsLocallyTrivial (Spec R) (tilde M).val := by
  refine isLocallyTrivial_of_pointwise _ (fun x => ?_)
  obtain ⟨g, hxg, ⟨e⟩⟩ := exists_away_linearEquiv (R : Type u) M x
  refine ⟨PrimeSpectrum.basicOpen g, hxg, ⟨?_⟩⟩
  letI : Module (Γ(Spec R, PrimeSpectrum.basicOpen g) : Type u)
      ((((restrictPresheafFunctor (Spec R) (PrimeSpectrum.basicOpen g)).obj
        (tilde M).val).obj (op (Over.mk (𝟙 (PrimeSpectrum.basicOpen g))))) : Type u) :=
    inferInstanceAs (Module (Γ(Spec R, PrimeSpectrum.basicOpen g) : Type u)
      (((tilde M).val.obj (op (PrimeSpectrum.basicOpen g))) : Type u))
  set s0 : ((((restrictPresheafFunctor (Spec R) (PrimeSpectrum.basicOpen g)).obj
      (tilde M).val).obj (op (Over.mk (𝟙 (PrimeSpectrum.basicOpen g))))) : Type u) :=
    generatorOf (sectionEquivScalar R M g e) with hs0
  refine (trivialIsoOfSectionPresieve (X := Spec R) (PrimeSpectrum.basicOpen g)
    ((restrictPresheafFunctor (Spec R) (PrimeSpectrum.basicOpen g)).obj (tilde M).val)
    (isSheaf_restrictModules _ (tilde M)) (isSheaf_unitOn _)
    (s := s0) (h := ?_)).symm
  intro W
  refine ⟨Presieve.functorPullback (Over.forget (PrimeSpectrum.basicOpen g))
      (fun (Z : (Spec R).Opens) (_ : Z ⟶ W.left) =>
        ∃ hh : (R : Type u), Z = PrimeSpectrum.basicOpen (hh * g)),
    overBasicPresieve_mem R g _ le_rfl W, ?_⟩
  rintro Z i ⟨hh, hZ⟩
  obtain ⟨Zl, Zr, Zhom⟩ := Z
  change Zl = PrimeSpectrum.basicOpen (hh * g) at hZ
  subst hZ
  have hgen : s0 = tildeAwayEquivScalar R M g (e.symm 1) := by
    rw [hs0]
    show (sectionEquivScalar R M g e).symm 1 = _
    show (tildeAwayEquivScalar R M g) ((awayEquivScalar R M g e).symm 1) = _
    congr 1
    show e.symm ((awayRingEquiv R g) 1) = _
    rw [map_one]
  have hkey2 : @restrictSec (Spec R) (PrimeSpectrum.basicOpen g)
      ((restrictPresheafFunctor (Spec R) (PrimeSpectrum.basicOpen g)).obj (tilde M).val) s0
      { left := PrimeSpectrum.basicOpen (hh * g), right := Zr, hom := Zhom }
      (overTermInst _ _)
      = tildeAwayEquivScalar R M (hh * g)
        (liftAwayMapA (R : Type u) g hh (M : Type u) (e.symm 1)) := by
    rw [hgen]
    exact DFunLike.congr_fun (tildeAwayEquiv_res R M g hh) (e.symm 1)
  rw [hkey2]
  exact bijective_smul_restrictGen R M g hh e

end

/-! ## ★出典の紐付け(`.src`) -/

def isLocallyTrivial_tilde.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——可逆加群の層は局所自明)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
