import ABC3.Found.GenEll.ComapMul
import ABC3.Found.GenEll.HeightConstruction

/-!
# [GenEll] Definition 1.2, (i) の足場 —— **引き戻しの底変換**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.5。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

## ★★★なぜこれが要るのか

原文の高さは **`ht_M̄ : X(ℚ̄) → ℝ`** である。
`X(ℚ̄)` の点は定義体 `F` を選んで `x_F : Spec 𝓞_F → X` で表すが、
★**`F` の取り方に依らない**ことを示さないと `ht` は well-defined でない。

★★`ArithDiv.lean` の `degNormalized_baseChange` が
「正規化次数は底変換で不変」を与えている。
★★★**足りないのは「引き戻したイデアルが底変換で拡大イデアルになる」ことである。**

## ★★これは `comap_mul` と同じ道具で出る

`ComapAffine.lean` の `ideal_comap_eq_map_of_isAffine`
(アフィンでは `comap` は `Ideal.map`)を `g : Spec 𝓞_K ⟶ Spec 𝓞_F` に当てるだけ。
★`Spec 𝓞_F`, `Spec 𝓞_K` はどちらもアフィンである。

★★★**`comap_mul` を取るために作った道具が、そのまま底変換に効いた。**
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField

variable (F K : Type) [Field F] [NumberField F] [Field K] [NumberField K]

variable {X : Scheme.{0}}

/-! ## ★★底変換の環準同型 -/

/-- ★**`g : Spec 𝓞_K ⟶ Spec 𝓞_F` が定める `𝓞_F → 𝓞_K`**。

★`ΓSpecIso` で両端を大域切断から環に戻す。 -/
noncomputable def baseRingHom (g : specRingOfIntegers K ⟶ specRingOfIntegers F) :
    CommRingCat.of (𝓞 F) ⟶ CommRingCat.of (𝓞 K) :=
  (Scheme.ΓSpecIso (CommRingCat.of (𝓞 F))).inv ≫ g.app ⊤ ≫
    (Scheme.ΓSpecIso (CommRingCat.of (𝓞 K))).hom

/-! ## ★★★引き戻しは底変換で拡大イデアルになる -/

set_option maxHeartbeats 2000000 in
/-- ★★★**引き戻したイデアルは、底を上げると拡大イデアルになる**。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

    `x_K^* D = (x_F^* D) · 𝓞_K`

★★これが `ht` を `X(ℚ̄)` の上で well-defined にするための有限素点側である。

★機構は `comap_comp` + `ideal_comap_eq_map_of_isAffine`(`ComapAffine.lean`)+
`Ideal.comap_symm`(同型に沿った `comap` は `map`)+ `Ideal.map_map`。
★★★**`comap_mul` を取るために作った道具がそのまま効いた。** -/
theorem pullbackIdeal_comp (D : X.IdealSheafData) (xF : specRingOfIntegers F ⟶ X)
    (g : specRingOfIntegers K ⟶ specRingOfIntegers F) :
    pullbackIdeal K D (g ≫ xF)
      = (pullbackIdeal F D xF).map (baseRingHom F K g).hom := by
  -- ★同型に沿った `comap` を `map` に直す(`.hom.hom` の形で述べる——
  --   `set` を使うと `↑e` になって後の照合が外れる)
  have keyF : ∀ J : Ideal Γ(specRingOfIntegers F, ⊤),
      Ideal.comap (Scheme.ΓSpecIso (CommRingCat.of (𝓞 F))).inv.hom J
        = J.map (Scheme.ΓSpecIso (CommRingCat.of (𝓞 F))).hom.hom :=
    fun J => Ideal.comap_symm
      (Scheme.ΓSpecIso (CommRingCat.of (𝓞 F))).commRingCatIsoToRingEquiv
  have keyK : ∀ J : Ideal Γ(specRingOfIntegers K, ⊤),
      Ideal.comap (Scheme.ΓSpecIso (CommRingCat.of (𝓞 K))).inv.hom J
        = J.map (Scheme.ΓSpecIso (CommRingCat.of (𝓞 K))).hom.hom :=
    fun J => Ideal.comap_symm
      (Scheme.ΓSpecIso (CommRingCat.of (𝓞 K))).commRingCatIsoToRingEquiv
  show Ideal.comap _ (Scheme.IdealSheafData.equivOfIsAffine (D.comap (g ≫ xF))) = _
  rw [Scheme.IdealSheafData.comap_comp]
  show Ideal.comap _ (((D.comap xF).comap g).ideal ⟨⊤, isAffineOpen_top _⟩) = _
  rw [ideal_comap_eq_map_of_isAffine, keyK]
  show _ = Ideal.map _ (Ideal.comap _
    (Scheme.IdealSheafData.equivOfIsAffine (D.comap xF)))
  rw [keyF]
  simp only [Ideal.map_map, Scheme.IdealSheafData.equivOfIsAffine_apply]
  -- ★`rw` は `⟨⊤, ⋯⟩` の証明項で照合が外れる。`Eq.trans` なら単一化で通る
  refine Eq.trans (Ideal.map_map _ _) ?_
  -- ★★規則: 圏の言葉のまま射の等式を作り、最後に 1 度だけ降りる
  have hmor : g.app ⊤ ≫ (Scheme.ΓSpecIso (CommRingCat.of (𝓞 K))).hom
      = (Scheme.ΓSpecIso (CommRingCat.of (𝓞 F))).hom ≫ baseRingHom F K g := by
    rw [baseRingHom, ← Category.assoc, ← Category.assoc, Iso.hom_inv_id, Category.id_comp]
  have hring := congrArg CommRingCat.Hom.hom hmor
  rw [CommRingCat.hom_comp, CommRingCat.hom_comp] at hring

  exact congrArg (fun φ => Ideal.map φ ((D.comap xF).ideal ⟨⊤, isAffineOpen_top _⟩)) hring

/-! ## ★出典の紐付け(`.src`)

★条つきである。`Definition 1.2, (i)` 全体には
アルキメデス側の底変換と `degNormalized_baseChange` との接続が要る。 -/

def pullbackIdeal_comp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Definition 1.2, (i)(引き戻しの底変換——有限素点側のみ)",
    sectionId := "genell-def-1-2-i" }

end ABC3.Found.GenEll
