import ABC3.Meta.Claim
import Mathlib.Algebra.Group.Equiv.Basic
import Mathlib.Algebra.Group.Int.Defs
import Mathlib.Algebra.Group.Units.Defs

/-!
# [Falt1] Theorem I.4.4 の posit(`Interface`)

原典: G. Faltings, *p-Adic Hodge Theory*, JAMS Vol.1 No.1 (1988) pp.255-299
(`papers.json` 短縮タグ `Falt1`)。**260 dpi 目視確認 2026-09-04**(物理 p.17-18、
印字 p.270-271、逐語の食い違い無し——数式・定理番号は完全に読める。地の文の
単語間スペース喪失のみ問題)。

`[LocProP] Lemma 2.1` が直接引用する箇所
(`Section I, Theorem 4.4, (i)`・`(iv)`)だけを posit する。

原文 (Falt1 p.271, (iv)):
> The morphism ΩR/V → H1(Δ, ρ-1R^(1)) given by the extension Eρ is
> induced by a functorial isomorphism
> H1(Δ,R^)/(p-torsion) ≅ ΩR/V(dlog∞) ⊗R (R⊗V V̄)^(-1),

★★**逸脱**: 原文は一般の `(R,{u_i})`(log 構造つきスキーム)について述べるが、
`[LocProP]` が消費するのは `H⁰`・`H¹ mod p-torsion` の**同型そのもの**だけなので、
`R`・`V`・`Δ` を抽象化し、対応する群の**名前だけ** posit する
(`H0`・`H1ModPTorsion`・`RVbar`・`OmegaTensor`)。定理の中身(同型の存在)は
posited data として直接受ける——証明していない。
-/

namespace ABC3.Interface.Falt1

universe u

/-- ★posit —— `[Falt1] Theorem I.4.4, (i)(iv)` が要る最小限の骨組み。 -/
structure PAdicHodgeSetup where
  /-- `Δ` = 関係する拡大の Galois 群。 -/
  Delta : Type u
  DeltaGrp : Group Delta
  /-- `H⁰(Δ, R̂(j))`。 -/
  H0 : ℤ → Type u
  H0Grp : ∀ j, AddCommGroup (H0 j)
  /-- `(R ⊗_V V̄)^(j)`。 -/
  RVbar : ℤ → Type u
  RVbarGrp : ∀ j, AddCommGroup (RVbar j)
  /-- **`[Falt1] Theorem I.4.4, (i)`**。 -/
  thm44i : ∀ j, H0 j ≃+ RVbar j
  /-- `H¹(Δ, R̂(j)) / (p-torsion)`(商そのものを posit する)。 -/
  H1ModPTorsion : ℤ → Type u
  H1ModPTorsionGrp : ∀ j, AddCommGroup (H1ModPTorsion j)
  /-- `Ω_{R/V}(dlog∞) ⊗_R (R⊗_V V̄)^(j)`。 -/
  OmegaTensor : ℤ → Type u
  OmegaTensorGrp : ∀ j, AddCommGroup (OmegaTensor j)
  /-- **`[Falt1] Theorem I.4.4, (iv)`**。 -/
  thm44iv : ∀ j, H1ModPTorsion j ≃+ OmegaTensor (j - 1)

/-- ★★**非退化性の witness**(具体項)—— すべて `ℤ` に潰す。 -/
@[reducible] def PAdicHodgeSetup.example : PAdicHodgeSetup.{0} where
  Delta := ℤˣ
  DeltaGrp := inferInstance
  H0 := fun _ => ℤ
  H0Grp := fun _ => inferInstance
  RVbar := fun _ => ℤ
  RVbarGrp := fun _ => inferInstance
  thm44i := fun _ => AddEquiv.refl ℤ
  H1ModPTorsion := fun _ => ℤ
  H1ModPTorsionGrp := fun _ => inferInstance
  OmegaTensor := fun _ => ℤ
  OmegaTensorGrp := fun _ => inferInstance
  thm44iv := fun _ => AddEquiv.refl ℤ

@[reducible] def PAdicHodgeSetup.nonvacuous : Nonempty (PAdicHodgeSetup.{0}) := ⟨PAdicHodgeSetup.example⟩

end ABC3.Interface.Falt1
