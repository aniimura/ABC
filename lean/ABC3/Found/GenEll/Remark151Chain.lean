/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.PullbackLocalization

/-!
# [GenEll] Remark 1.5.1 —— **鎖をつなぐ**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

## ★★★★★★★★★★ここまでの 6 本を 1 本に繋ぐ

| ファイル | 出したもの |
|---|---|
| `IsoDescent.lean` | 生成ファイバーの同型が有限段へ降りる |
| `DivisorDescent.lean` | 因子ごと降りる（可換正方形） |
| `ConductorDescent.lean` | 正方形 ⟹ 段 `n` で因子が一致（`comap_eq_of_square`） |
| `LocalizationBridge.lean` | `𝓞_F[1/N]` で一致 ⟹ `N` を割らない `v` で導手係数が一致 |
| `PullbackLocalization.lean` | 引き戻しイデアルを任意の環の上で |
| ★**本ファイル** | ★**点の持ち上げを与えれば導手が一致する** |

## ★★★★★★「`Σ` の外で」の正体

原文の「`Σ` の外の素点 `v` で」は、形式化では**2 つの仮定**に分かれる:

1. `hM : ∀ m ∈ M, m ∉ v.asIdeal` —— `v` が `M`（＝ `N` の冪）を避ける。
2. `xA : Spec 𝓞_F[1/N] ⟶ X_n` —— **点が段 `n` のモデルへ持ち上がる**。

★★2 は `𝓞_F[1/N]` が `ℤ[1/n!]`-代数であることから作れる（`N = n!` と取る）が、
`ratTower n = Subring.closure {(n!)⁻¹}` に `IsLocalization` の characterization が
まだ無いので、**本ファイルでは仮定として受ける**。
★★★受け取り側（`Remark 1.5.1` の本体）にとっては、
「そういう `v` に対して成り立つ」という形なので**影響しない**。
-/

namespace ABC3.Found.GenEll

open CategoryTheory Limits AlgebraicGeometry

variable (F : Type) [Field F] [NumberField F]

/-- ★★**任意の環の上でも「点を動かすか因子を動かすかは同じ」**。 -/
theorem pullbackIdealOf_target (B : CommRingCat.{0}) {X X' : Scheme.{0}}
    (D' : X'.IdealSheafData) (φ : X ⟶ X') (x : Spec B ⟶ X) :
    pullbackIdealOf B D' (x ≫ φ) = pullbackIdealOf B (D'.comap φ) x := by
  simp [pullbackIdealOf, Scheme.IdealSheafData.comap_comp]

/-- ★★★★★★★★★★**[GenEll] Remark 1.5.1 —— 鎖の全体**。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

`DivisorDescent.lean` が出す「有限段 `n` での同型＋可換正方形」と、
**その段へ点が持ち上がること**（`xA`）から、
`Σ` の外の素点 `v` で**導手の係数が一致する**。

★★★これが `Skeleton/GenEll/Section1.lean` の `remark_1_5_1` が受けている
仮定 `hagree` そのものである。

★機構は 3 本の合成:
`comap_eq_of_square`（幾何）→ `pullbackIdealOf_target`（点と因子の入れ替え）
→ `conductorADiv_fin_eq_of_localized`（局所化の橋）。 -/
theorem conductorADiv_fin_eq_of_lift
    (M : Submonoid (NumberField.RingOfIntegers F)) (A : Type) [CommRing A]
    [Algebra (NumberField.RingOfIntegers F) A] [IsLocalization M A]
    (v : FinitePlace F) (hM : ∀ m ∈ M, m ∉ v.asIdeal)
    {Z Z' X X' : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (f' : X' ⟶ Spec (CommRingCat.of ℤ))
    (iZ : Z ⟶ X) (iZ' : Z' ⟶ X') [IsClosedImmersion iZ] [IsClosedImmersion iZ']
    {n : ℕᵒᵖ}
    (φ : bcObj f n ⟶ bcObj f' n) (ψ : bcObj (iZ ≫ f) n ⟶ bcObj (iZ' ≫ f') n)
    [IsIso φ] [IsIso ψ]
    (hsq : ψ ≫ bcBC (iZ' ≫ f') f' iZ' n = bcBC (iZ ≫ f) f iZ n ≫ φ)
    (xA : Spec (CommRingCat.of A) ⟶ bcObj f n)
    (xF : specRingOfIntegers F ⟶ X) (xF' : specRingOfIntegers F ⟶ X')
    (hlift : Spec.map (CommRingCat.ofHom (algebraMap (NumberField.RingOfIntegers F) A)) ≫ xF
      = xA ≫ pullback.snd (overRatTowerDiagram.obj n).hom f)
    (hlift' : Spec.map (CommRingCat.ofHom (algebraMap (NumberField.RingOfIntegers F) A)) ≫ xF'
      = xA ≫ φ ≫ pullback.snd (overRatTowerDiagram.obj n).hom f')
    (hI : pullbackIdeal F iZ.ker xF ≠ 0) (hJ : pullbackIdeal F iZ'.ker xF' ≠ 0) :
    (conductorADiv F iZ.ker xF).fin v = (conductorADiv F iZ'.ker xF').fin v := by
  refine conductorADiv_fin_eq_of_localized F M A v hM iZ.ker iZ'.ker xF xF' hI hJ ?_
  rw [hlift, hlift', pullbackIdealOf_target, ← Category.assoc, pullbackIdealOf_target,
    comap_eq_of_square f f' iZ iZ' φ ψ hsq]
  exact (pullbackIdealOf_target (CommRingCat.of A) _ φ xA).symm

/-! ### ★出典の紐付け(`.src`)

★★**項目全体の `.src` はまだ置かない。** 残っているのは
「`𝓞_F[1/n!]` が `ratTower n`-代数であること」から `xA` を実際に作る段だけである
——`ratTower n = Subring.closure {(n!)⁻¹}` に `IsLocalization` の
characterization がまだ無い。 -/

def conductorADiv_fin_eq_of_lift.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(hagree——点の持ち上げを仮定として受ける形)",
    sectionId := "genell-rem-1-5-1" }

def conductorADiv_fin_eq_of_lift.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "comap_eq_of_square(正方形 ⟹ 段 n で因子が一致)"
      (.inProject "ABC3" "ABC3.Found.GenEll.comap_eq_of_square") 9,
    .citation "[ABC3]" "conductorADiv_fin_eq_of_localized(局所化の橋)"
      (.inProject "ABC3" "ABC3.Found.GenEll.conductorADiv_fin_eq_of_localized") 9,
    .citation "[ABC3]" "exists_pair_descent(対 (X, D) の spreading out——正方形を出す)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_pair_descent") 9,
    .implicitStep
      ("★★原文の「Σ の外の素点 v で」は形式化では 2 つの仮定に分かれる: " ++
       "(1) v が M(＝N の冪)を避ける、(2) 点が段 n のモデルへ持ち上がる(xA)。" ++
       "★(2) は 𝓞_F[1/N] が ℤ[1/n!]-代数であることから作れるが、" ++
       "ratTower n = Subring.closure {(n!)⁻¹} に IsLocalization の " ++
       "characterization がまだ無いので本ファイルでは仮定として受けた") 9,
    .implicitStep
      ("★★★残る段: ratTower n が ℤ の powers (n!) による局所化であることを示し、" ++
       "その普遍性で ratTower n ⟶ 𝓞_F[1/n!] を作る。" ++
       "★そこから pullback.lift で xA が出る") 9 ]

end ABC3.Found.GenEll
