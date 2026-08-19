/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.ProL.LPart
import ABC3.Found.ProL.FinitePrimary

/-!
# `M ≅ ∏_l M[l]` —— 可換副有限群の pro-`l` 分解(チェーン `prol` の葉 `decomposition`)

原文 (FrdI p.52):
> Thus, M decomposes as a direct product of pro-l groups M [l], where l varies over

## ★★★測って選んだ道(2026-08-19)

★**逆向き `Φ : ∏_l M[l] → M` を先に作る**。
`Φ y` は「各有限商で `∏_{p ∣ |M/U|} y_p` を取る」整合族の一意な持ち上げである
(`existsUnique_mk_eq_of_compatible`)。

★★**連続性は圏論を経由しない** —— 副有限群では
「開正規部分群の剰余類」が開基をなす(`exist_openNormalSubgroup_sub_open_nhds_of_one`)ので、
`y ↦ f y U`(**有限**積、離散空間への写像)が連続であることだけで足りる。

★★★全単射を示したら、**コンパクト→Hausdorff の連続全単射は同相**
(`Continuous.homeoOfEquivCompactToT2`)で終わる。

| 段 | 中身 |
|---|---|
| `lComp` | `x` の `l`-成分(各商で `primProj`、整合性は `primProj_naturality`) |
| `Φ` | `∏_l M[l] → M` の構成と群準同型性 |
| 単射 | `primProj` を当てて成分ごとに分離 |
| 全射 | `lComp` で作った族が `x` に戻る(`prod_primProj`) |
| 連続 | 剰余類が開基／`f y U` が有限積 |
-/

namespace ABC3.Found.ProL

open CategoryTheory

universe u

/-- ★`ℕ` から `Nat.Primes` への全域写像(素数でなければ `2`)。 -/
noncomputable def toPrimes (p : ℕ) : Nat.Primes :=
  if h : p.Prime then ⟨p, h⟩ else ⟨2, Nat.prime_two⟩

@[simp] theorem toPrimes_val {p : ℕ} (h : p.Prime) : (toPrimes p).1 = p := by
  simp [toPrimes, h]

section CommProfinite

variable {M : Type u} [CommGroup M] [TopologicalSpace M] [IsTopologicalGroup M]
  [CompactSpace M] [TotallyDisconnectedSpace M]

/-- ★`M` を `ProfiniteGrp` として見る。 -/
abbrev asProfiniteGrp (M : Type u) [CommGroup M] [TopologicalSpace M] [IsTopologicalGroup M]
    [CompactSpace M] [TotallyDisconnectedSpace M] : ProfiniteGrp.{u} := ProfiniteGrp.of M

/-! ## ★1. `l`-成分 -/

variable (M) in
/-- ★`x` の `l`-成分を与える整合族。 -/
noncomputable def lCompFam (l : Nat.Primes) (x : M) (U : OpenNormalSubgroup M) :
    M ⧸ U.toSubgroup :=
  primProj (M ⧸ U.toSubgroup) l.1 (QuotientGroup.mk x)

variable (M) in
theorem lCompFam_compatible (l : Nat.Primes) (x : M) :
    ∀ (U V : OpenNormalSubgroup M) (h : U.toSubgroup ≤ V.toSubgroup),
      QuotientGroup.map U.toSubgroup V.toSubgroup (MonoidHom.id M) h (lCompFam M l x U)
        = lCompFam M l x V := by
  intro U V h
  show QuotientGroup.map _ _ _ h (primProj _ l.1 (QuotientGroup.mk x))
    = primProj _ l.1 (QuotientGroup.mk x)
  rw [primProj_naturality _ l.2]
  congr 1

variable (M) in
/-- ★★**`x` の `l`-成分**。 -/
noncomputable def lComp (l : Nat.Primes) (x : M) : M :=
  (existsUnique_mk_eq_of_compatible (asProfiniteGrp M) (lCompFam M l x)
    (lCompFam_compatible M l x)).choose

variable (M) in
theorem lComp_spec (l : Nat.Primes) (x : M) (U : OpenNormalSubgroup M) :
    (QuotientGroup.mk (lComp M l x) : M ⧸ U.toSubgroup)
      = primProj (M ⧸ U.toSubgroup) l.1 (QuotientGroup.mk x) :=
  (existsUnique_mk_eq_of_compatible (asProfiniteGrp M) (lCompFam M l x)
    (lCompFam_compatible M l x)).choose_spec.1 U

variable (M) in
/-- ★`l`-成分は `M[l]` に入る。 -/
theorem lComp_mem (l : Nat.Primes) (x : M) : lComp M l x ∈ lPart M l.1 := by
  intro U
  have h := lComp_spec M l x U
  obtain ⟨k, hk⟩ := primProj_mem_primaryComponent (A := M ⧸ U.toSubgroup) l.2
    (QuotientGroup.mk x)
  exact ⟨k, by rw [h]; exact hk⟩

end CommProfinite

end ABC3.Found.ProL
