import Mathlib.FieldTheory.Finite.GaloisField
import Mathlib.FieldTheory.PrimitiveElement

/-!
# 有限体上には任意の次数のモニック既約多項式がある

論文にも我々のモデルにも固有でない、**一般の**結果。mathlib へ出せる形で書く
(`Found/ResidueFieldFinite.lean` と同じ位置づけ)。

## なぜ要るのか

`Found/PGC/UnramifiedExtension.lean` の
`exists_isUnramifiedAdjoin_of_irreducible` は「剰余体上のモニック既約
多項式 `g`」から次数 `deg g` の**不分岐**拡大を作る。したがって
「各次数 `n` に不分岐拡大が存在する」を言うには、有限体上に次数 `n` の
モニック既約多項式があることだけが足りない。

## 証明の筋(mathlib の道具だけで数行)

`F` の標数を `p`、`m := [F : ZMod p]` とする。`GaloisField p (m*n)` は
`[· : ZMod p] = m*n` の有限体で、`m ∣ m*n` だから
`FiniteField.nonempty_algHom_of_finrank_dvd` が `F` からの `ZMod p`-代数
準同型を与える。これで `F`-代数とみなすと、塔の公式
(`Module.finrank_mul_finrank`)から `[GaloisField p (m*n) : F] = n`。
有限体は完全体なので分離拡大、`Field.exists_primitive_element` で
原始元 `α` が取れ、`minpoly F α` が求めるモニック既約 `n` 次多項式。
-/

namespace ABC3.Found

open Polynomial

/-- **有限体上には任意の次数 `n ≥ 1` のモニック既約多項式がある**。 -/
theorem exists_monic_irreducible_natDegree_eq (F : Type*) [Field F] [Finite F] (n : ℕ)
    (hn : n ≠ 0) :
    ∃ g : Polynomial F, g.Monic ∧ Irreducible g ∧ g.natDegree = n := by
  haveI := Fintype.ofFinite F
  obtain ⟨p, hp⟩ := CharP.exists F
  haveI : Fact (Nat.Prime p) := ⟨CharP.char_is_prime F p⟩
  letI := ZMod.algebra F p
  have hm : Module.finrank (ZMod p) F ≠ 0 := by
    have h1 := FiniteField.pow_finrank_eq_card p F
    intro h
    rw [h, pow_zero] at h1
    have := Fintype.one_lt_card (α := F)
    omega
  have hL : Module.finrank (ZMod p) (GaloisField p (Module.finrank (ZMod p) F * n))
      = Module.finrank (ZMod p) F * n :=
    GaloisField.finrank p (by simpa using ⟨hm, hn⟩)
  obtain ⟨φ⟩ := FiniteField.nonempty_algHom_of_finrank_dvd
    (F := ZMod p) (K := F) (L := GaloisField p (Module.finrank (ZMod p) F * n))
    (by rw [hL]; exact Dvd.intro n rfl)
  letI : Algebra F (GaloisField p (Module.finrank (ZMod p) F * n)) := φ.toRingHom.toAlgebra
  haveI : IsScalarTower (ZMod p) F (GaloisField p (Module.finrank (ZMod p) F * n)) :=
    IsScalarTower.of_algebraMap_eq (fun z => (φ.commutes z).symm)
  haveI : FiniteDimensional F (GaloisField p (Module.finrank (ZMod p) F * n)) :=
    Module.Finite.of_finite
  have htower := Module.finrank_mul_finrank (ZMod p) F
    (GaloisField p (Module.finrank (ZMod p) F * n))
  rw [hL] at htower
  have hrk : Module.finrank F (GaloisField p (Module.finrank (ZMod p) F * n)) = n :=
    Nat.eq_of_mul_eq_mul_left (Nat.pos_of_ne_zero hm) htower
  obtain ⟨α, hα⟩ := Field.exists_primitive_element F
    (GaloisField p (Module.finrank (ZMod p) F * n))
  have hint : IsIntegral F α := IsIntegral.of_finite F α
  refine ⟨minpoly F α, minpoly.monic hint, minpoly.irreducible hint, ?_⟩
  rw [← IntermediateField.adjoin.finrank hint, hα, IntermediateField.finrank_top', hrk]


/-- **次数が一致する有限体拡大には既約多項式の根がある**——`AdjoinRoot g`
(次数 `deg g` の有限体)から `E` への `F`-代数準同型が
`FiniteField.nonempty_algHom_of_finrank_dvd` で取れる。 -/
theorem exists_root_of_finrank_eq (F : Type*) [Field F] [Finite F] (g : Polynomial F)
    (hgi : Irreducible g) (E : Type*) [Field E] [Finite E] [Algebra F E]
    (hdeg : Module.finrank F E = g.natDegree) :
    ∃ b : E, Polynomial.aeval b g = 0 := by
  haveI : Fact (Irreducible g) := ⟨hgi⟩
  have h1 : Module.finrank F (AdjoinRoot g) = g.natDegree := by
    rw [(AdjoinRoot.powerBasis (f := g) hgi.ne_zero).finrank, AdjoinRoot.powerBasis_dim]
  obtain ⟨φ⟩ := FiniteField.nonempty_algHom_of_finrank_dvd (F := F) (K := AdjoinRoot g) (L := E)
    (by rw [h1, hdeg])
  refine ⟨φ (AdjoinRoot.root g), ?_⟩
  rw [Polynomial.aeval_algHom_apply, AdjoinRoot.aeval_eq, AdjoinRoot.mk_self, map_zero]

end ABC3.Found
