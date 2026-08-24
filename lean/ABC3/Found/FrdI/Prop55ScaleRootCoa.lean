/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop55ScaleRootInv
import ABC3.Found.FrdI.Prop55CongrIso
import ABC3.Found.FrdI.Cor410Birat

/-!
# [FrdI] Proposition 5.5, (ii) —— `Σ_k` は co-angular pre-step を保つ

原文 (FrdI p.104):
> (ii) There is a natural equivalence of categories [compatible with the functors to the

★`Prop55ScaleRootInv.lean` で `Σ_k` が `IsLinear` / `IsBaseIsomorphism` /
`IsPreStep` / `IsIsometric` を**保存し反射する**ことを閉じた。
本ファイルはそれを在庫の `isCoAngular_map_of_equiv`(`Theorem 3.4, (iii)`)へ流して
**`coaPreProp` の保存**(`psiBiratCor` の入力 `hfwd`)を得る。

## ★★道筋

`isCoAngular_map_of_equiv` が要求するのは**擬逆 `e.inverse` についての 4 型の保存**である。
★これは counit の同型四角形

```
Σ_k (e.inverse.map g) ≫ counit.app Y = counit.app X ≫ g
```

を `Prop55CongrIso.lean` の `*_congr_iso` で渡り、
`Prop55ScaleRootInv.lean` の**反射**(`*_of_scaleRootHom`)で降ろせばよい。

★★★これで `psiBiratCor`(`Corollary 4.10`)に `Σ_k` を流せる ——
`Proposition 5.5, (ii)` の「一般の根を根 1 に落とす」4 段のうち
**3 段目(`Σ_k` を birat へ降ろす)の入力が揃った**。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section ScaleRootCoa

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}

/-- ★★★★★**`Σ_k` は co-angular 性を保つ**。

★在庫の `isCoAngular_map_of_equiv`(`Theorem 3.4, (iii)`)に、
counit の四角形を渡る `*_congr_iso` と `Σ_k` の**反射**を組み合わせる。 -/
theorem isCoAngular_scaleRootHom (k : ℕ+) {X Y : PfRootObj P F} (f : HomRoot P F X Y)
    (h : IsCoAngular (pfRootPre P F) f) :
    IsCoAngular (pfRootPre P F) (scaleRootHom (F := F) k f) := by
  set e := (scaleRootFunctor P F k).asEquivalence with he
  have hsq : ∀ {X' Y' : PfRootObj P F} (g : X' ⟶ Y'),
      e.functor.map (e.inverse.map g) ≫ (e.counitIso.app Y').hom
        = (e.counitIso.app X').hom ≫ g := fun {_ _} g => e.counitIso.hom.naturality g
  refine isCoAngular_map_of_equiv e ?_ ?_ ?_ ?_ f h
  · intro X' Y' g hg
    exact isLinear_of_scaleRootHom k
      ((isLinear_congr_iso (pfRootPre P F) _ g
        (e.counitIso.app X') (e.counitIso.app Y') (hsq g)).mpr hg)
  · intro X' Y' g hg
    exact isIsometric_of_scaleRootHom k
      ((isIsometric_congr_iso (pfRootPre P F) _ g
        (e.counitIso.app X') (e.counitIso.app Y') (hsq g)).mpr hg)
  · intro X' Y' g hg
    exact isPreStep_of_scaleRootHom k
      ((isPreStep_congr_iso (pfRootPre P F) _ g
        (e.counitIso.app X') (e.counitIso.app Y') (hsq g)).mpr hg)
  · intro X' Y' g hg
    exact isBaseIsomorphism_of_scaleRootHom k
      ((isBaseIsomorphism_congr_iso (pfRootPre P F) _ g
        (e.counitIso.app X') (e.counitIso.app Y') (hsq g)).mpr hg)

/-- ★★★★★★**`psiBiratCor` の入力 `hfwd`** —— `Σ_k` は co-angular pre-step を保つ。

★これで `Corollary 4.10` の `psiBiratCor` に `Σ_k` を流せる。 -/
theorem coaPreProp_scaleRootHom (k : ℕ+) {X Y : PfRootObj P F} (f : HomRoot P F X Y)
    (h : coaPreProp (pfRootPre P F) f) :
    coaPreProp (pfRootPre P F) ((scaleRootFunctor P F k).map f) :=
  ⟨isCoAngular_scaleRootHom k f h.1, isPreStep_scaleRootHom k h.2⟩

/-- ★★★★★**`Σ_k` は co-angular 性を反射する** —— 在庫 `isCoAngular_of_map` に
`Σ_k` の**保存**(`Prop55ScaleRootInv.lean`)を渡すだけ。 -/
theorem isCoAngular_of_scaleRootHom (k : ℕ+) {X Y : PfRootObj P F} (f : HomRoot P F X Y)
    (h : IsCoAngular (pfRootPre P F) (scaleRootHom (F := F) k f)) :
    IsCoAngular (pfRootPre P F) f :=
  isCoAngular_of_map (scaleRootFunctor P F k)
    (fun _ hg => isLinear_scaleRootHom k hg)
    (fun _ hg => isIsometric_scaleRootHom k hg)
    (fun _ hg => isPreStep_scaleRootHom k hg)
    (fun _ hg => isBaseIsomorphism_scaleRootHom k hg)
    f h

/-- ★★★★★★**`psiBiratCor` の入力 `hbwd`** —— `Σ_k` は co-angular pre-step を反射する。 -/
theorem coaPreProp_of_scaleRootHom (k : ℕ+) {X Y : PfRootObj P F} (f : HomRoot P F X Y)
    (h : coaPreProp (pfRootPre P F) ((scaleRootFunctor P F k).map f)) :
    coaPreProp (pfRootPre P F) f :=
  ⟨isCoAngular_of_scaleRootHom k f h.1, isPreStep_of_scaleRootHom k h.2⟩

/-! ## ★2. `Σ_k` を birat 化へ降ろす -/

/-- ★★★★★★**[FrdI] Proposition 5.5, (ii) の 3 段目** ——
**根の一斉倍化 `Σ_k` を `(𝒞^pf)^birat` へ降ろした関手**。

★在庫 `psiBiratCor`(`Corollary 4.10`)に `Σ_k` と `hfwd` を流すだけ。 -/
noncomputable def scaleRootBirat (k : ℕ+) (Gpf : Frobenioid (pfRootPre P F)) :
    BiratCat (pfRootPre P F) Gpf ⥤ BiratCat (pfRootPre P F) Gpf :=
  psiBiratCor Gpf Gpf (scaleRootFunctor P F k)
    (fun {_ _} f hf => coaPreProp_scaleRootHom k f hf)

/-- ★★★★★★**降ろした `Σ_k` も圏同値**。

★★これが `Proposition 5.5, (ii)` の「一般の根を根 1 に落とす」段である。 -/
theorem scaleRootBirat_isEquivalence (k : ℕ+) (Gpf : Frobenioid (pfRootPre P F)) :
    (scaleRootBirat (F := F) k Gpf).IsEquivalence :=
  psiBiratCor_isEquivalence Gpf Gpf (scaleRootFunctor P F k)
    (fun {_ _} f hf => coaPreProp_scaleRootHom k f hf)
    (fun {_ _} f hf => coaPreProp_of_scaleRootHom k f hf)

/-- ★★★★★★**`Σ_k` の birat 版が与える Hom の全単射**。

★`Proposition 5.5, (ii)` の 4 段の 3 段目そのもの ——
`Hom(⟨A,r⟩,⟨B,s⟩) ≃ Hom(⟨A,k·r⟩,⟨B,k·s⟩)`。 -/
noncomputable def scaleRootBiratHomEquiv (k : ℕ+) (Gpf : Frobenioid (pfRootPre P F))
    (X Y : BiratCat (pfRootPre P F) Gpf) :
    (X ⟶ Y) ≃ ((scaleRootBirat (F := F) k Gpf).obj X
      ⟶ (scaleRootBirat (F := F) k Gpf).obj Y) :=
  haveI := scaleRootBirat_isEquivalence (F := F) k Gpf
  Equiv.ofBijective (fun f : X ⟶ Y => (scaleRootBirat (F := F) k Gpf).map f)
    ⟨fun _ _ h => (scaleRootBirat (F := F) k Gpf).map_injective h,
      fun g => (scaleRootBirat (F := F) k Gpf).map_surjective g⟩

end ScaleRootCoa

/-! ### ★出典の紐付け -/

/-- ★★★★★★locator —— `Proposition 5.5, (ii)` の梱包に要る
「根の一斉倍化 `Σ_k` は co-angular pre-step を保つ」(`psiBiratCor` の入力)。 -/
def coaPreProp_scaleRootHom.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — Σ_k は co-angular pre-step を保つ",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
