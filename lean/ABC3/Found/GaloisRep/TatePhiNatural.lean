/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.NormRepNatural
import ABC3.Found.GaloisRep.TatePhi

/-!
# Tate 一意化の**係数の自然性**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★`galois-equivariant-tate-uniformization` の段 2b

台帳の 7 段のうち:

| 段 | 内容 | 状態 |
|---|---|---|
| 1 | 固定した完備体での**選ばれた**同型 | ★済（`tatePhiAddEquivAll`） |
| 2a | `normRep` の自然性 | ★済（`NormRepNatural.lean`） |
| **2b** | **`tateAOf` / `tateWOf` の自然性** | ★**本ファイル** |
| 2c | `tatePtPair` の自然性（級数と `σ` の交換） | 未 |
| 3 | 有限拡大に対する自然性 | 未 |
| 4 | `K̄` への帰納極限 | 未 |
| 5 | `l`-捩れの完全列 | ★消費側は済 |

## ★★★★★★機構は 2 本

`tateAOf S c` は `.choose` だが、**`tateAOf_spec` が特徴づける**:

    algebraMap R K (tateAOf S c) = normRep S.v S.Q S.hQ c

★段 2a（`normRep` の自然性）と `S.hinj`（`algebraMap R K` の単射性、
`TateSetup` の欄）を繋ぐと `tateAOf` の自然性が出る。

★★`tateWOf` は `tateAOf · tateWOf = q` から**消去**で出る
（`tateAOf ≠ 0` は `normRep` が単元だから）。

★★★`σR S.q = S.q` は仮説ではなく**導出される**——
`algebraMap S.q = S.Q` と `σU S.Q = S.Q` と単射性から。

## ★残っている段（明示）

★`tatePhi` 全体の同変性には `tatePtPair` の自然性（段 2c）が要る。
★★そこは**級数と `σ` を交換させる段**であり、`σ` の連続性
（付値の拡大の一意性から等距離写像）が前提である。
-/

namespace ABC3.Found.GaloisRep

variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [Algebra R K]

/-! ## ★係数は零でない -/

/-- ★**`tateAOf` は零でない**（`normRep` が単元だから）。 -/
theorem tateAOf_ne_zero (S : TateSetup R I K) (c : Kˣ ⧸ Subgroup.zpowers S.Q) :
    tateAOf S c ≠ 0 := by
  intro h
  have h1 := (tateAOf_spec S c).1
  rw [h, map_zero] at h1
  exact (normRep S.v S.Q S.hQ c).ne_zero h1.symm

/-! ## ★★母数は固定される -/

/-- ★★**`σR S.q = S.q` は導出される**（仮説ではない）。

★`algebraMap R K S.q = S.Q` と `σU S.Q = S.Q` と `S.hinj` から。 -/
theorem tateSetup_q_map (S : TateSetup R I K)
    (σR : R →+* R) (σK : K →+* K)
    (hcompat : ∀ r : R, σK (algebraMap R K r) = algebraMap R K (σR r))
    (σU : Kˣ →* Kˣ) (hσU : ∀ u : Kˣ, ((σU u : Kˣ) : K) = σK (u : K))
    (hσq : σU S.Q = S.Q) : σR S.q = S.q := by
  refine S.hinj ?_
  rw [← hcompat, S.hQq, ← hσU, hσq]

/-! ## ★★★★★★★係数の自然性 -/

/-- ★★★★★★★**`tateAOf` は `σ` で自然**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

> `σR (tateAOf S c) = tateAOf S (σ_* c)`

★★`tateAOf` は `.choose` だが `tateAOf_spec` が特徴づけるので、
段 2a（`normRep` の自然性）と単射性を繋ぐだけで出る。 -/
theorem tateAOf_map (S : TateSetup R I K)
    (σR : R →+* R) (σK : K →+* K)
    (hcompat : ∀ r : R, σK (algebraMap R K r) = algebraMap R K (σR r))
    (σU : Kˣ →* Kˣ) (hσU : ∀ u : Kˣ, ((σU u : Kˣ) : K) = σK (u : K))
    (hσq : σU S.Q = S.Q) (hσv : ∀ u, vAdd S.v (σU u) = vAdd S.v u)
    (c : Kˣ ⧸ Subgroup.zpowers S.Q) :
    σR (tateAOf S c)
      = tateAOf S (QuotientGroup.map _ _ σU (zpowers_le_comap_self S.Q σU hσq) c) := by
  refine S.hinj ?_
  rw [← hcompat, (tateAOf_spec S c).1, ← hσU,
    normRep_map_self S.v S.Q S.hQ σU hσq hσv c]
  exact ((tateAOf_spec S _).1).symm

/-- ★★★★★★★**`tateWOf` も `σ` で自然**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★`tateAOf · tateWOf = q` の両辺に `σR` を当て、`tateAOf` の自然性で
左因子を揃えてから**消去**する（`R` は整域、`tateAOf ≠ 0`）。 -/
theorem tateWOf_map (S : TateSetup R I K)
    (σR : R →+* R) (σK : K →+* K)
    (hcompat : ∀ r : R, σK (algebraMap R K r) = algebraMap R K (σR r))
    (σU : Kˣ →* Kˣ) (hσU : ∀ u : Kˣ, ((σU u : Kˣ) : K) = σK (u : K))
    (hσq : σU S.Q = S.Q) (hσv : ∀ u, vAdd S.v (σU u) = vAdd S.v u)
    (c : Kˣ ⧸ Subgroup.zpowers S.Q) :
    σR (tateWOf S c)
      = tateWOf S (QuotientGroup.map _ _ σU (zpowers_le_comap_self S.Q σU hσq) c) := by
  set c' := QuotientGroup.map (Subgroup.zpowers S.Q) (Subgroup.zpowers S.Q) σU
    (zpowers_le_comap_self S.Q σU hσq) c with hc'
  have hprod : ∀ d : Kˣ ⧸ Subgroup.zpowers S.Q, tateAOf S d * tateWOf S d = S.q := by
    intro d
    refine S.hinj ?_
    rw [(tateAOf_spec S d).2.2, S.hQq]
  have hA : σR (tateAOf S c) = tateAOf S c' := tateAOf_map S σR σK hcompat σU hσU hσq hσv c
  have h1 : σR (tateAOf S c) * σR (tateWOf S c) = S.q := by
    rw [← map_mul, hprod c]
    exact tateSetup_q_map S σR σK hcompat σU hσU hσq
  rw [hA] at h1
  have h2 : tateAOf S c' * tateWOf S c' = S.q := hprod c'
  have hne : tateAOf S c' ≠ 0 := tateAOf_ne_zero S c'
  exact mul_left_cancel₀ hne (h1.trans h2.symm)

/-! ### ★出典の紐付け(`.src`)

★★`Definition 3.3`（Tate 一意化）の**同変性の段 2b** であって、
`tatePhi` 全体の同変性ではない。 -/

def tateAOf_map.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——係数 tateAOf の自然性)",
    sectionId := "genell-def-3-3" }

def tateWOf_map.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——係数 tateWOf の自然性)",
    sectionId := "genell-def-3-3" }

def tateAOf_map.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "normRep_map_self(正規化代表の同変性——段 2a)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.normRep_map_self") 15,
    .citation "[ABC3]" "tateAOf_spec(tateAOf を normRep で特徴づける)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateAOf_spec") 15,
    .implicitStep
      ("★★★.choose で定義されたものでも、**特徴づけがあれば自然性は出る**。" ++
       "★tateAOf は algebraMap R K (tateAOf S c) = normRep … で決まり、" ++
       "algebraMap R K は単射(TateSetup の hinj)である") 15,
    .implicitStep
      ("★★残る段 2c: tatePtPair の自然性(級数と σ の交換)。" ++
       "★前提は σ の連続性であり、付値の拡大の一意性から等距離写像になる") 15 ]

end ABC3.Found.GaloisRep
