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

## ★★★★★★用途によっては関手に束ねなくてよい(2026-08-25 の測定)

★★**`Proposition 5.5, (iii)` の `𝒞^pf` が birationally Frobenius-normalized** を出すには、
在庫の `isFrobeniusNormalized_transport`(`TypeTransport.lean`)を使う。
★これが要求するのは **1 つの対象での `End` の乗法同型 ＋ 3 つの対応**
(`IsBaseIdentity` / `degFr` / `OTri`)だけであり、**関手は要らない**。

★★★下の `biratPfHomEquivRoot` を `A = B`、`n = m`、`k = n·n` で読むと

    End_{(𝒞^birat)^pf}(biratUp (rtObj A n))  ≃  End_{(𝒞^pf)^birat} ⟨A, n⟩

がそのまま出る(2026-08-25 に型検査で確認)。
★したがって下の「難所」——`Hom(X,Y)` と `Hom(Y,Z)` で共役の同型が別物になる問題——は
**単一対象の `End` では起きない**(共役の同型が 1 つに固定されるから)。

★残るのは (1) この全単射が**乗法的**であること
(`biratPfHom_id` / `biratPfHom_comp` ＋ `scaleRootBirat` が関手 ＋ 同型による共役)、
(2) 3 つの対応(`biratBase_biratPfHom` / `biratDeg_biratPfHom` / `otimes_biratPfHom` ＋
`rootBase_scaleRootHom` / `rootDeg_scaleRootHom`)の 2 つだけである。

## ★★★★関手に束ねるときの難所(2026-08-25 の測定、次に着手する人へ)

★この全単射を**関手**にするには `map_id` / `map_comp` が要るが、
第 3 段の共役に使う同型

```
⟨A,n⟩ ≅ ⟨A^{(m)}, k⟩   (k = n·m)
```

は **`k` を通して「相手の対象の根 `m`」に依存する**。
★したがって `Hom(X,Y)` と `Hom(Y,Z)` で使う同型が**別物**になり、
共役をそのまま合成しても `Hom(X,Z)` の共役にならない。

★★**逃げ道 A(直接)**: 在庫の `compRoot_eq_lift`
(「**合成は共通倍数の根で計算してよい**」、`Prop55BiratPf.lean` の ★ に記述)で
3 対象を共通の根 `n·m·l` に揃えてから比べる。
`scaleRootFunctor` の `map_comp` がまさにこの形で書かれているので、同じ骨が使える。

★★**逃げ道 B(普遍性、2026-08-25 に見つけた別ルート)**:
**向きを逆にして**、`𝒞^birat` の**普遍性**で作る。

```
Ω : 𝒞^pf ⥤ (𝒞^birat)^pf     (代表元を toBiratCat で送るだけ。★在庫 idxToBirat がある)
```

は co-angular pre-step を**同型に送る**(birat では同型だから)。
★もし「co-angular pre-step を同型に送る関手は `𝒞^birat` を経由する」という
**birat 化の普遍性**があれば、`Θ : (𝒞^pf)^birat ⥤ (𝒞^birat)^pf` が
**関手として無料で出る**(合成の coherence を手で書かずに済む)。

★★ただし在庫の `HomBirat` は**添字圏の余極限**として作ってあり、
`HomBirat.mk` / `exists_rep` / `sound` はあるが **`HomBirat.desc`(普遍性)は無い**。
★どちらの逃げ道も同程度の作業量に見えるが、**B は再利用が効く**
(`Corollary 4.10` の `psiBiratCor` を手で組んでいる所も普遍性で書き直せる)。

## ★★★★★訂正(2026-08-25)—— 根の「等式」は**取れない**

前の版では「`F'` の根を `F` と揃える(`rtObj (biratPre P G) F' A d = rtObj P F A d`)」と
書いたが、★**これは仮定としても成り立たない**。
`rtObj Q F' A d = (F'.frobDegSurj A d).choose` であり、`Exists` は `Prop` なので
`choose` は**命題だけで決まる** —— `𝒞` の中の命題と `𝒞^birat` の中の命題は別物だから、
どちらの `choose` も動かせない。

★★**正しい手当ては `rootIso`(根の不変性)である**。
`rtObj (biratPre P G) F' A d` と `toBiratCat` で送った `rtObj P F A d` は
**どちらも `A` の次数 `d` の Frobenius 拡大**なので、
`F'.frobDegUniq` が(構造射と両立する)同型を与える。
その同型に沿って `rootIso` を当てれば `HomPf` が移る ——
★等式は要らず、**同型で十分**である。

★★★これは `Corollary 5.4` の継ぎ目で `modelType_equiv` の
`Classical.choice` を `thm_5_2_iv` の明示的な切断に置き換えたのと**同じ型の話**だが、
結論は逆向きである ——
そこでは切断を**揃えられた**(`BaseSection.map` があった)のに対し、
ここは揃えられないので**同型で逃げる**。 -/
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

/-! ## ★4. `F'` の根と `F` の根を **`rootIso` で**繋ぐ -/

/-- ★★★★**`F'` の根から `F` の根への同型** ——
`rtObj (biratPre P G) F' A d` と `toBiratCat` で送った `rtObj P F A d` は
どちらも `A` の次数 `d` の Frobenius 拡大なので、`F'.frobDegUniq` が同型を与える。

★★等式は取れない(`Exists.choose` は命題だけで決まる)が、**同型で十分**である。 -/
theorem exists_rtObj_birat_iso {G : Frobenioid P} (F' : FrobenioidCore (biratPre P G))
    (Z : C) (d : ℕ+) :
    ∃ β : rtObj (biratPre P G) F' (biratUp P G Z) d ⟶ biratUp P G (rtObj P F Z d),
      IsIso β := by
  have hfrob : IsFrobeniusType (biratPre P G) ((toBiratCat P G).map (rtExt P F Z d)) :=
    (birat_isFrobeniusType_iff P G (rtExt P F Z d)).mpr
      ⟨(rtExt_frobType P F Z d).1.1, (rtExt_frobType P F Z d).2⟩
  have hdeg : (biratPre P G).degFr (rtExt (biratPre P G) F' (biratUp P G Z) d)
      = (biratPre P G).degFr ((toBiratCat P G).map (rtExt P F Z d)) := by
    rw [rtExt_degFr]
    show (d : ℕ+) = biratDeg (toHomBirat (P := P) (G := G) (rtExt P F Z d))
    rw [biratDeg_toHomBirat, rtExt_degFr]
  obtain ⟨β, hβ, -⟩ := F'.frobDegUniq _ _ _
    (rtExt (biratPre P G) F' (biratUp P G Z) d)
    ((toBiratCat P G).map (rtExt P F Z d))
    (rtExt_frobType (biratPre P G) F' (biratUp P G Z) d) hfrob hdeg
  exact ⟨β, hβ⟩

/-- ★★★★★★★**[FrdI] Proposition 5.5, (ii) の射の全単射(一般の根、仮定なし)** ——

`Hom_{(𝒞^birat)^pf}((A,n),(B,m)) ≃ Hom_{(𝒞^pf)^birat}(⟨A,n⟩,⟨B,m⟩)`。

★`F'` の根と `F` の根の差は `rootIso`(根の不変性)で吸収する ——
**等式は要らない**。 -/
noncomputable def biratPfHomEquivRootFull (hfi : IsOfFrobeniusIsotropicType P)
    (hiso : ∀ X : C, IsIsotropic P X) {G : Frobenioid P}
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G))
    (A B : C) (n m k : ℕ+) (hk : k * 1 = n * m) (hk' : k * 1 = m * n) :
    HomRoot (biratPre P G) F'
        (⟨biratUp P G A, n⟩ : PfRootObj (biratPre P G) F') ⟨biratUp P G B, m⟩
      ≃ HomBirat (pfRootPre P F) Gpf (⟨A, n⟩ : PfRootObj P F) ⟨B, m⟩ :=
  haveI hA := (exists_rtObj_birat_iso (F := F) F' A m).choose_spec
  haveI hB := (exists_rtObj_birat_iso (F := F) F' B n).choose_spec
  ((rootIso (F := F') (exists_rtObj_birat_iso (F := F) F' A m).choose
      (isFrobeniusType_of_isIso (biratPre P G) _)
      (exists_rtObj_birat_iso (F := F) F' B n).choose
      (isFrobeniusType_of_isIso (biratPre P G) _)
      (by rw [isLinear_of_isIso (biratPre P G) _, isLinear_of_isIso (biratPre P G) _])).symm
    |>.toEquiv).trans
    (biratPfHomEquivRoot hfi hiso Gpf F' A B n m k hk hk')

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

/-- ★★★★★★★locator —— `Proposition 5.5, (ii)` の
**一般の根での射の全単射(仮定なし版)** —— `F'` の根と `F` の根の差は
`rootIso`(根の不変性)で吸収する。 -/
def biratPfHomEquivRootFull.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — 一般の根での射の全単射(F' の根は rootIso で吸収)",
    sectionId := "frdi-prop-5-5" }

/-- ★★★★★★locator —— `Corollary 4.10` の `Ψ^birat` に `Σ_k` を流したもの
(`Proposition 5.5, (ii)` の 3 段目)。 -/
def scaleRootBirat.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — Σ_k を (𝒞^pf)^birat へ降ろす(圏同値)",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
