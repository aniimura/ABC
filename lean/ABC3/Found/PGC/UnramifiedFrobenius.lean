import ABC3.Found.PGC.UnramifiedRootsOfUnity

/-!
# `Gal(K(x)/K)` の整数環・剰余体への作用(`adjoinField` 側)

目標は「Frobenius は `p` と素な 1 の冪根に `q` 乗で作用する」
——`q` を群論的に読む(LCFT 非経由)道の核心の算術——だが、
本ファイルではその**手前まで**しか進んでいない。何が壁かを記録する。

## できたこと

`σ : Gal(K(x)/K)` を **`𝒪[(adjoinField K x).carrier]` の環同型**として持ち上げ
(`integersEquivOf`)、値の計算(`integersEquivOf_val`)まで通した。
さらに剰余体の同型(`residueEquivOf`)も作れる。

★`algEquivIntegers`(`adjoinIntegers` 側)との合成で作ろうとすると
**kernel が止まる**。`algEquivIntegers` と同じ手順で**直接**書くと通る
(`tools/lean-idioms.md` #59)。

## ★★★測定: `adjoinIntegers` と `𝒪[(adjoinField K x).carrier]` の境界は越えられない

`exists_frobenius`(`UnramifiedExtension.lean`)が与える Frobenius の性質は
**`ResidueField (adjoinIntegers K x)` 上**で述べられている。それを
`𝓀[(adjoinField K x).carrier]` 側に運ぶには

```
(integersEquivAdjoinIntegers K x w).1 = w.1
```

が要る。これは定義上ほとんど自明(`toFun z := ⟨z.1, _⟩`)だが、
両辺の**型**が `↥K⟮x⟯` と `(adjoinField K x).carrier` で異なるため、
`rfl` も `Subtype.ext rfl` も

* 既定の heartbeats: `(kernel) deterministic timeout`(212 秒)
* `maxHeartbeats 1000000`: `(deterministic) timeout at whnf`(126 秒)

で**通らない**(2026-09-05 実測)。`𝒪[·]` が `Valuation.integer` +
スペクトルノルム由来の `Valued` インスタンスで、`whnf` が展開しきれない。

**したがって、この橋は「書き方を工夫する」では越えられない。**
越えるには `exists_frobenius` 自体を `𝒪[(adjoinField K x).carrier]` 側で
述べ直す(`UnramifiedExtension.lean` の該当部分を書き換える)しかない。
次にここへ戻るときは、その書き換えから始めること。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped Valued

variable {p : ℕ} [Fact p.Prime]

/-- `σ` を `𝒪[(adjoinField K x).carrier]` の環同型として。

★`algEquivIntegers`(`adjoinIntegers` 側)との合成で作ると
中間体 2 層をまたぐ `rfl` が kernel を止める(`tools/lean-idioms.md` #59)。
`algEquivIntegers` と同じ手順で**直接**書く。 -/
noncomputable def integersEquivOf (K : PAdicLocalField p) (x : K.closure)
    (σ : (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) :
    𝒪[(adjoinField K x).carrier] ≃+* 𝒪[(adjoinField K x).carrier] where
  toFun z := ⟨σ (z : (adjoinField K x).carrier),
    (mem_adjoinIntegers_iff_mem_integers K x (σ (z : (adjoinField K x).carrier))).mpr
      (by
        show ‖σ (z : (adjoinField K x).carrier)‖ ≤ 1
        rw [norm_algEquiv K x σ]
        exact (mem_adjoinIntegers_iff_mem_integers K x (z : (adjoinField K x).carrier)).mp z.2)⟩
  invFun z := ⟨σ.symm (z : (adjoinField K x).carrier),
    (mem_adjoinIntegers_iff_mem_integers K x (σ.symm (z : (adjoinField K x).carrier))).mpr
      (by
        show ‖σ.symm (z : (adjoinField K x).carrier)‖ ≤ 1
        rw [norm_algEquiv K x σ.symm]
        exact (mem_adjoinIntegers_iff_mem_integers K x (z : (adjoinField K x).carrier)).mp z.2)⟩
  left_inv z := by apply Subtype.ext; simp
  right_inv z := by apply Subtype.ext; simp
  map_mul' a b := by
    apply Subtype.ext
    show σ ((a : (adjoinField K x).carrier) * (b : (adjoinField K x).carrier))
      = σ (a : (adjoinField K x).carrier) * σ (b : (adjoinField K x).carrier)
    exact map_mul σ _ _
  map_add' a b := by
    apply Subtype.ext
    show σ ((a : (adjoinField K x).carrier) + (b : (adjoinField K x).carrier))
      = σ (a : (adjoinField K x).carrier) + σ (b : (adjoinField K x).carrier)
    exact map_add σ _ _

@[simp] theorem integersEquivOf_val (K : PAdicLocalField p) (x : K.closure)
    (σ : (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
    (w : 𝒪[(adjoinField K x).carrier]) :
    ((integersEquivOf K x σ w : 𝒪[(adjoinField K x).carrier]) :
      (adjoinField K x).carrier) = σ (w : (adjoinField K x).carrier) := rfl

/-- `σ` が誘導する剰余体 `𝓀[(adjoinField K x).carrier]` の環同型。 -/
noncomputable def residueEquivOf (K : PAdicLocalField p) (x : K.closure)
    (σ : (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) :
    𝓀[(adjoinField K x).carrier] ≃+* 𝓀[(adjoinField K x).carrier] :=
  IsLocalRing.ResidueField.mapEquiv (integersEquivOf K x σ)

end ABC3.Found.PGC
