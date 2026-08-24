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
  [IsConnected D]

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

/-! ## ★3. 4 段を繋ぐ —— **一般の根での射の全単射** -/

/-- ★★★★★★★**[FrdI] Proposition 5.5, (ii) の射の全単射(一般の根)** ——

```
Hom_{(𝒞^birat)^pf}((A,n),(B,m)) = Hom^pf_{𝒞^birat}(A^{(m)}, B^{(n)})   -- 定義
    ≃ Hom_{(𝒞^pf)^birat}(⟨A^{(m)},1⟩, ⟨B^{(n)},1⟩)   -- biratPfHomEquiv
    ≃ Hom_{(𝒞^pf)^birat}(⟨A^{(m)},k⟩, ⟨B^{(n)},k⟩)   -- scaleRootBiratHomEquiv
    ≃ Hom_{(𝒞^pf)^birat}(⟨A,n⟩, ⟨B,m⟩)               -- pfRoot_exists_iso_root
```

★★これが `Prop55BiratPfFun.lean` が「根 1 の場合しか無い」と書いていた所の一般化である。
★左辺の `A^{(m)}` は **`F` の根**で書いてある —— `F'` の根と揃える手当ては
呼び出し側(`rtObj (biratPre P G) F' = rtObj P F`)で行う。

## ★★★★関手に束ねるときの難所(2026-08-25 の測定、次に着手する人へ)

★この全単射を**関手**にするには `map_id` / `map_comp` が要るが、
第 3 段の共役に使う同型

```
⟨A,n⟩ ≅ ⟨A^{(m)}, k⟩   (k = n·m)
```

は **`k` を通して「相手の対象の根 `m`」に依存する**。
★したがって `Hom(X,Y)` と `Hom(Y,Z)` で使う同型が**別物**になり、
共役をそのまま合成しても `Hom(X,Z)` の共役にならない。

★★**逃げ道**: 在庫の `compRoot_eq_lift`
(「**合成は共通倍数の根で計算してよい**」、`Prop55BiratPf.lean` の ★ に記述)で
3 対象を共通の根 `n·m·l` に揃えてから比べる。
`scaleRootFunctor` の `map_comp` がまさにこの形で書かれているので、同じ骨が使える。

★もう 1 つの手当ては `F'` の根を `F` と揃えること
(`rtObj (biratPre P G) F' A d = rtObj P F A d`)——
★★これは `Corollary 5.4` の継ぎ目で**基準切断を揃えた**のと同じ型の手当てである
(`Cor54SeamCls.lean` / `Cor54SeamSec.lean` を見よ)。 -/
noncomputable def biratPfHomEquivRoot (hfi : IsOfFrobeniusIsotropicType P)
    (hiso : ∀ X : C, IsIsotropic P X) {G : Frobenioid P}
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G))
    (A B : C) (n m k : ℕ+) (hk : k * 1 = n * m) (hk' : k * 1 = m * n) :
    HomPf (biratPre P G) F' (biratUp P G (rtObj P F A m)) (biratUp P G (rtObj P F B n))
      ≃ HomBirat (pfRootPre P F) Gpf (⟨A, n⟩ : PfRootObj P F) ⟨B, m⟩ :=
  haveI := (pfRoot_exists_iso_root (F := F) A n m (k * 1) hk).choose_spec
  haveI := (pfRoot_exists_iso_root (F := F) B m n (k * 1) hk').choose_spec
  (biratPfHomEquiv hfi hiso Gpf F' (rtObj P F A m) (rtObj P F B n)).trans
    ((scaleRootBiratHomEquiv (F := F) k Gpf
        (⟨rtObj P F A m, 1⟩ : PfRootObj P F) (⟨rtObj P F B n, 1⟩ : PfRootObj P F)).trans
      (Iso.homCongr
        ((toBiratCat (pfRootPre P F) Gpf).mapIso
          (asIso (pfRoot_exists_iso_root (F := F) A n m (k * 1) hk).choose)).symm
        ((toBiratCat (pfRootPre P F) Gpf).mapIso
          (asIso (pfRoot_exists_iso_root (F := F) B m n (k * 1) hk').choose)).symm))

end ScaleRootCoa

/-! ### ★出典の紐付け -/

/-- ★★★★★★locator —— `Proposition 5.5, (ii)` の梱包に要る
「根の一斉倍化 `Σ_k` は co-angular pre-step を保つ」(`psiBiratCor` の入力)。 -/
def coaPreProp_scaleRootHom.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — Σ_k は co-angular pre-step を保つ",
    sectionId := "frdi-prop-5-5" }

/-- ★★★★★★★locator —— `Proposition 5.5, (ii)` の**一般の根での射の全単射**
(4 段の鎖を繋いだもの)。 -/
def biratPfHomEquivRoot.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — 一般の根での射の全単射(4 段の鎖)",
    sectionId := "frdi-prop-5-5" }

/-- ★★★★★★locator —— `Corollary 4.10` の `Ψ^birat` に `Σ_k` を流したもの
(`Proposition 5.5, (ii)` の 3 段目)。 -/
def scaleRootBirat.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — Σ_k を (𝒞^pf)^birat へ降ろす(圏同値)",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
