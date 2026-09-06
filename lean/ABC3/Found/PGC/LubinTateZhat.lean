import ABC3.Found.PGC.AbelianDecomposition
import ABC3.Found.PGC.UnramifiedZhat

/-!
# 経路 Λ4 ∘ Λ5 —— `Gal(K_π·K^ur/K) ≅ 𝒪_K^× × Ẑ`

[pGC] Proposition 1.1 の壁 `TorsionCyclotomeIsCyclotomic` に向けた、
局所類体論の**アーベル部分**の分解を 1 本にまとめる。

## なぜ別ファイルか

`Found/PGC/UnramifiedZhat.lean`(第 1047)は `Ẑ` との同型を
**剰余体を経由せずに**出しており、Lubin–Tate 一式に依存していない。
そこへ `AbelianDecomposition` を import すると、
★**Lubin–Tate 塔(12,367 行)が `Ẑ` の葉に乗ってしまう**。
合成だけを欲しい消費者のために、依存を足す側を別ファイルに切り出す。

## 材料(いずれも 2026-09-06、すべて `sorry` 0)

* `AbelianDecomposition.lean` の `exists_lubinTateUnramified_decomposition`
  —— `K_π ⊓ K^ur = ⊥` と `Gal(K_π·K^ur/K) ≃* 𝒪_K^× × Gal(K^ur/K)`
* `UnramifiedZhat.lean` の `unramifiedClosureGalEquivZHat`
  —— `Gal(K^ur/K) ≃ₜ* Ẑ`

## ☆逸脱の記録

無し。両者の合成のみで、新しい仮定も読み替えも入れていない。
★`Ẑ` は `Multiplicative ℤ` の副有限完備化として与えている
(mathlib の `ProfiniteCompletion` が `GrpCat` の上にあるため)。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

variable {p : ℕ} [Fact p.Prime]

/-- ★★**`Gal(K_π·K^ur/K) ≅ 𝒪_K^× × Ẑ`**。

経路 Λ の Λ4(アーベル分解)と Λ5(`Gal(K^ur/K) ≅ Ẑ`)の合成である。

★部分群・`π`・`f` はすべて `∃` の内側に閉じ込めてある
(結論の statement に自由なパラメータを出さない —— 2026-09-06 の 10 例目の退化を避ける形)。 -/
theorem exists_lubinTateUnramified_decomposition_zhat (K : PAdicLocalField p) :
    ∃ E : IntermediateField K.carrier K.closure,
      Normal K.carrier E ∧ E ⊓ unramifiedClosure K = ⊥ ∧
      Nonempty ((E ≃ₐ[K.carrier] E) ≃* (𝒪[K.carrier])ˣ) ∧
      Nonempty (((E ⊔ unramifiedClosure K : IntermediateField K.carrier K.closure) ≃ₐ[K.carrier]
          (E ⊔ unramifiedClosure K : IntermediateField K.carrier K.closure)) ≃*
        ((𝒪[K.carrier])ˣ × ZHat)) := by
  obtain ⟨E, hn, hinf, h1, ⟨e⟩⟩ := exists_lubinTateUnramified_decomposition K
  exact ⟨E, hn, hinf, h1,
    ⟨e.trans ((MulEquiv.refl _).prodCongr (unramifiedClosureGalEquivZHat K).toMulEquiv)⟩⟩

end ABC3.Found.PGC
