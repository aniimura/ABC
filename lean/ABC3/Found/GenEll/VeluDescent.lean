/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.VeluGaloisInvariant
import Mathlib.FieldTheory.Galois.Infinite
import Mathlib.FieldTheory.IsAlgClosed.AlgebraicClosure
import ABC3.Meta.Claim

/-!
# 第 1154 ブロック —— **`Gal`-安定な部分群による Vélu の商は `L` に降りる**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★これは何か——`LCyclicReading` の節点 1 の**完成**

第 1153 で「座標の集合 `S` が `Φ` で保たれるなら Vélu の商の係数は `Φ` で固定される」
（`fixesCoeffs_veluQuotientFull`）が取れた。★残っていた 1 歩は

    `L̄^{Gal(L̄/L)} = L`  すなわち  「すべての `σ` で固定 ⟹ `L` の元」

である。☆mathlib に**そのものがあった**——`InfiniteGalois.mem_range_algebraMap_iff_fixed`
（`Mathlib/FieldTheory/Galois/Infinite.lean`）。

    `x ∈ Set.range (algebraMap k K) ↔ ∀ f : Gal(K/k), f x = x`   （`[IsGalois k K]`）

★標数 0 なら `IsGalois L (AlgebraicClosure L)` はインスタンスで出る
（分離性は `CharZero`、正規性は代数閉包）。

## ★★★★★★★★到達点

    `veluQuotientFull_descends`

——`W : WeierstrassCurve L` と、`L̄` の座標の有限集合 `S` が**すべての `σ` で保たれる**なら、
`veluQuotientFull (W⁄L̄) S` は **`L` 上の曲線の底変換**である。

★★★これで `Skeleton/GenEll/LCyclicReading.lean` の**節点 1 が閉じた**。
☆残るのは節点 2（`Lemma 3.5` を安定直線の側で述べ直す）と節点 3（`Lemma 3.7` の言い直し）。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep WeierstrassCurve

/-! ## ★★★★★★★★すべての `σ` で固定なら `L` の元 -/

variable {L : Type} [Field L] [CharZero L]

local notation "Lbar" => AlgebraicClosure L

/-- ★★★★★★★★**`L̄` の元がすべての `σ` で固定されるなら `L` の元**——★**無条件**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆mathlib の `InfiniteGalois.mem_range_algebraMap_iff_fixed` そのものである。
★標数 0 なので `IsGalois L L̄` はインスタンスで出る。 -/
theorem mem_range_of_fixed (x : Lbar)
    (h : ∀ σ : Lbar ≃ₐ[L] Lbar, σ x = x) :
    x ∈ Set.range (algebraMap L Lbar) :=
  (InfiniteGalois.mem_range_algebraMap_iff_fixed x).2 h

/-! ## ★★★★★★★★★★★★底変換した曲線の係数は `σ` で固定される -/

/-- ☆`L` 上の曲線を `L̄` へ上げたものの係数は、どの `σ` でも固定される。 -/
theorem fixesCoeffs_baseChange (W : WeierstrassCurve L) (σ : Lbar ≃ₐ[L] Lbar) :
    FixesCoeffs (W.map (algebraMap L Lbar)).toAffine (σ : Lbar ≃+* Lbar) :=
  ⟨σ.commutes W.a₁, σ.commutes W.a₂, σ.commutes W.a₃, σ.commutes W.a₄, σ.commutes W.a₆⟩

/-! ## ★★★★★★★★★★★★★★★★★★★★商は `L` に降りる -/

/-- ★★★★★★★★★★★★★★★★★★★★
**`Gal`-安定な座標集合による Vélu の商は `L` 上の曲線の底変換**——★**無条件**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★これが `Skeleton/GenEll/LCyclicReading.lean` の**節点 1 そのもの**である。
☆原文の `H_L ⊆ E_L` は `l`-巡回部分群**スキーム**であり、その点は `L` 有理とは限らない。
★だが `H` が `Gal`-安定なら `E/H` は `L` 上で定義される——それが本定理である。

☆`S` に「部分群の座標集合であること」は要求していない。必要なのは
**`Gal` の作用で保たれること**だけである。 -/
theorem veluQuotientFull_descends (W : WeierstrassCurve L)
    (S : Finset (Lbar × Lbar))
    (hS : ∀ σ : Lbar ≃ₐ[L] Lbar, ∀ z ∈ S, semiPair (σ : Lbar ≃+* Lbar) z ∈ S) :
    ∃ W' : WeierstrassCurve L,
      W'.map (algebraMap L Lbar)
        = veluQuotientFull (W.map (algebraMap L Lbar)) S := by
  set WB : WeierstrassCurve Lbar := W.map (algebraMap L Lbar) with hWB
  -- ★商の `a₄`・`a₆` はすべての `σ` で固定される
  have hfix : ∀ σ : Lbar ≃ₐ[L] Lbar,
      FixesCoeffs (veluQuotientFull WB S).toAffine (σ : Lbar ≃+* Lbar) := fun σ =>
    fixesCoeffs_veluQuotientFull WB (fixesCoeffs_baseChange W σ) S (hS σ)
  obtain ⟨a₄', ha₄'⟩ : (veluQuotientFull WB S).a₄ ∈ Set.range (algebraMap L Lbar) :=
    mem_range_of_fixed _ (fun σ => (hfix σ).a₄)
  obtain ⟨a₆', ha₆'⟩ : (veluQuotientFull WB S).a₆ ∈ Set.range (algebraMap L Lbar) :=
    mem_range_of_fixed _ (fun σ => (hfix σ).a₆)
  refine ⟨⟨W.a₁, W.a₂, W.a₃, a₄', a₆'⟩, ?_⟩
  have h₁ : (veluQuotientFull WB S).a₁ = algebraMap L Lbar W.a₁ := rfl
  have h₂ : (veluQuotientFull WB S).a₂ = algebraMap L Lbar W.a₂ := rfl
  have h₃ : (veluQuotientFull WB S).a₃ = algebraMap L Lbar W.a₃ := rfl
  refine WeierstrassCurve.ext ?_ ?_ ?_ ?_ ?_
  · exact h₁.symm
  · exact h₂.symm
  · exact h₃.symm
  · exact ha₄'
  · exact ha₆'

/-! ## ★★★★★★★★★★★★★★★★降下はさらなる底変換で保たれる -/

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★
**降下した商はどの底変換でも Vélu の商のまま**——★**無条件**（第 1179）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`W'⁄L̄ = veluQuotientFull (W⁄L̄) S`（第 1154）を `f : L̄ →+* M` で送ると、
`W'⁄M = veluQuotientFull (W⁄M) (f(S))` になる。

★★★これが `Skeleton/GenEll/LCyclicReading.lean` の**節点 2d の中身**である
——`M` を `L̄_v` と取れば `E' ⊗ L_v` を局所の Vélu の商と合わせる `hcurveEq` が出る。
☆`veluQuotientFull_map`（第 679、在庫）をそのまま使う。 -/
theorem veluQuotientFull_descends_map (W W' : WeierstrassCurve L)
    (S : Finset (Lbar × Lbar))
    (hW' : W'.map (algebraMap L Lbar) = veluQuotientFull (W.map (algebraMap L Lbar)) S)
    {M : Type} [Field M] (f : Lbar →+* M) :
    W'.map (f.comp (algebraMap L Lbar))
      = veluQuotientFull (W.map (f.comp (algebraMap L Lbar)))
          (S.image (fun Q => (f Q.1, f Q.2))) := by
  rw [← WeierstrassCurve.map_map, hW', veluQuotientFull_map, ← WeierstrassCurve.map_map]

/-! ## ★★★★★★★★★★★★★★★★★★★★任意の代数閉体へ運ぶ -/

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★
**`Gal`-安定な商はどの代数閉体の上でも Vélu の商**——★**無条件**（第 1180）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`M` が `L` を含む代数閉体なら `IsAlgClosed.lift` で `L̄ ↪ M` が取れるので、
第 1154（降下）と第 1179（底変換）を合わせて

    `W'⁄M = veluQuotientFull (W⁄M) (ι(S))`

が出る。★★★`M` を `L̄_v`（`L_v` の代数閉包）と取れば、
これが `LCyclicReading` の節点 2d が要求する `hcurveEq` の形である。
☆残るのは `ι(S)` が Tate 一意化の下で `μ_l` の座標集合になること
（`Lemma 3.2, (i)`）を当てる段だけである。 -/
theorem veluQuotientFull_descends_algClosed (W : WeierstrassCurve L)
    (S : Finset (Lbar × Lbar))
    (hS : ∀ σ : Lbar ≃ₐ[L] Lbar, ∀ z ∈ S, semiPair (σ : Lbar ≃+* Lbar) z ∈ S)
    (M : Type) [Field M] [IsAlgClosed M] [Algebra L M] :
    ∃ (W' : WeierstrassCurve L) (f : Lbar →+* M),
      (∀ x : L, f (algebraMap L Lbar x) = algebraMap L M x) ∧
      W'.map (algebraMap L M)
        = veluQuotientFull (W.map (algebraMap L M))
            (S.image (fun Q => (f Q.1, f Q.2))) := by
  obtain ⟨W', hW'⟩ := veluQuotientFull_descends W S hS
  letI ι : Lbar →ₐ[L] M := IsAlgClosed.lift
  refine ⟨W', ι.toRingHom, fun x => ι.commutes x, ?_⟩
  have hcomp : (ι.toRingHom : Lbar →+* M).comp (algebraMap L Lbar) = algebraMap L M := by
    ext x
    exact ι.commutes x
  have h := veluQuotientFull_descends_map W W' S hW' (ι.toRingHom : Lbar →+* M)
  rwa [hcomp] at h

/-! ## ★★★★★★★★★★★★`L(H)`——座標の定義体を取る -/

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★
**有限個の座標は有限次拡大の中に入る**——★**無条件**（第 1195）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`S` の座標で生成する中間体 `L(S)` は `L` 上有限次である
（`L̄` は `L` 上代数的だから）。

★★★これが `Skeleton/GenEll/LCyclicReading.lean` の節点 2d-1 の
**`L'' ≔ L(H)` を取る段**である（第 1194 の道の第 1 段）。 -/
theorem exists_finite_subextension (S : Finset (Lbar × Lbar)) :
    ∃ M : IntermediateField L Lbar, FiniteDimensional L M ∧
      ∀ z ∈ S, z.1 ∈ M ∧ z.2 ∈ M := by
  classical
  set T : Set Lbar := (↑(S.image Prod.fst) : Set Lbar) ∪ (↑(S.image Prod.snd) : Set Lbar) with hT
  have hTfin : T.Finite := (S.image Prod.fst).finite_toSet.union (S.image Prod.snd).finite_toSet
  haveI : Finite T := hTfin
  have hint : ∀ x ∈ T, IsIntegral L x := fun x _ => Algebra.IsIntegral.isIntegral x
  refine ⟨IntermediateField.adjoin L T, IntermediateField.finiteDimensional_adjoin hint, ?_⟩
  intro z hz
  refine ⟨IntermediateField.subset_adjoin L T (Or.inl ?_),
    IntermediateField.subset_adjoin L T (Or.inr ?_)⟩
  · exact Finset.mem_coe.2 (Finset.mem_image.2 ⟨z, hz, rfl⟩)
  · exact Finset.mem_coe.2 (Finset.mem_image.2 ⟨z, hz, rfl⟩)

/-! ## ★★★★★★★★曲線の底変換は単射 -/

/-- ★★★★★★★★**単射な環準同型での底変換は曲線の等式を落とす**——★**無条件**（第 1196）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`W.map f = W'.map f` かつ `f` が単射なら `W = W'`。
★★これが `L̄` での等式を中間体 `M` へ**降ろす**ための一行であり、
第 1194 の道の**段 2**（`E' ⊗ L'' = veluQuotientFull (E ⊗ L'') …`）で使う。 -/
theorem weierstrassCurve_map_injective {R A : Type*} [CommRing R] [CommRing A] (f : R →+* A)
    (hf : Function.Injective f) {W W' : WeierstrassCurve R} (h : W.map f = W'.map f) :
    W = W' := by
  refine WeierstrassCurve.ext ?_ ?_ ?_ ?_ ?_
  · exact hf (congrArg WeierstrassCurve.a₁ h)
  · exact hf (congrArg WeierstrassCurve.a₂ h)
  · exact hf (congrArg WeierstrassCurve.a₃ h)
  · exact hf (congrArg WeierstrassCurve.a₄ h)
  · exact hf (congrArg WeierstrassCurve.a₆ h)

/-! ## ★出典の紐付け(`.src`) -/

def mem_range_of_fixed.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(L̄ の元がすべての σ で固定されるなら L の元。★無条件)",
    sectionId := "genell-lemma-3-5" }

def mem_range_of_fixed.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "InfiniteGalois.mem_range_algebraMap_iff_fixed(L̄^Gal = L)"
      (.inMathlib "InfiniteGalois.mem_range_algebraMap_iff_fixed") 1 ]

def fixesCoeffs_baseChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(底変換した曲線の係数はどの σ でも固定される。★無条件)",
    sectionId := "genell-lemma-3-5" }

def weierstrassCurve_map_injective.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(単射な環準同型での底変換は曲線の等式を落とす。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_finite_subextension.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(有限個の座標は有限次拡大の中に入る——L(H) を取る段。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluQuotientFull_descends_algClosed.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Gal-安定な商はどの代数閉体の上でも Vélu の商。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluQuotientFull_descends_algClosed.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "IsAlgClosed.lift(代数的な体は代数閉体に埋め込める)"
      (.inMathlib "IsAlgClosed.lift") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1180）**——`M` を `L̄_v`（`L_v` の代数閉包）と取れば、" ++
       "これが `LCyclicReading` の節点 2d が要求する `hcurveEq` の形である。" ++
       "☆残るのは `ι(S)` が Tate 一意化の下で `μ_l` の座標集合になること" ++
       "（`Lemma 3.2, (i)`）を当てる段だけである。") 1 ]

def veluQuotientFull_descends_map.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(降下した商はどの底変換でも Vélu の商のまま。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluQuotientFull_descends_map.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "veluQuotientFull_map(第 679、在庫)"
      (.inProject "ABC3" "ABC3.Found.GenEll.veluQuotientFull_map") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1179）**——`LCyclicReading` の**節点 2d の中身**である。" ++
       "☆`M` を `L̄_v` と取れば、`E' ⊗ L_v` を局所の Vélu の商と合わせる " ++
       "`hcurveEq` が出る。★残るのは `L̄ ↪ L̄_v` を `IsAlgClosed.lift` で取り、" ++
       "安定直線の像が `μ_l` になること（`Lemma 3.2, (i)`）を当てる段だけである。") 2 ]

def veluQuotientFull_descends.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Gal-安定な座標集合による Vélu の商は L 上の曲線の底変換。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluQuotientFull_descends.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "fixesCoeffs_veluQuotientFull(第 1153、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.fixesCoeffs_veluQuotientFull") 1,
    .citation "[mathlib]" "InfiniteGalois.mem_range_algebraMap_iff_fixed(L̄^Gal = L)"
      (.inMathlib "InfiniteGalois.mem_range_algebraMap_iff_fixed") 1,
    .implicitStep
      ("★★★★**2026-09-01（第 1154）**——`Skeleton/GenEll/LCyclicReading.lean` の" ++
       "**節点 1 が閉じた**。☆原文の `H_L` は部分群スキームでその点は `L` 有理とは限らないが、" ++
       "`Gal`-安定なら `E/H` は `L` 上で定義される。" ++
       "★残るのは節点 2（`Lemma 3.5` を安定直線の側で述べ直す）と" ++
       "節点 3（`Lemma 3.7` の言い直し）である。") 4 ]

end ABC3.Found.GenEll
