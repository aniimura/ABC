/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.Thm38KummerExists
import ABC3.Found.GenEll.Thm38Alpha
import ABC3.Found.GaloisRep.TateGaloisStab
import ABC3.Found.GenEll.EllModuliGalois
import ABC3.Meta.Claim

/-!
# 第 1160 ブロック —— **`α` が像に入る段の橋**（`Skeleton`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19–p.20。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★これは何か

`Theorem 3.8` を `Found` で閉じるのに要るものは**ちょうど 2 つ**である。

| # | 中身 | 状態 |
|---|---|---|
| I | `l`-巡回の読みの橋（`HasLCyclicVelu` ⟷ `HasLCyclicJ`） | ☆`Skeleton/GenEll/LCyclicReading.lean`（進捗枠 **2.7 / 3**、残り重み 4） |
| II | **`α` が mod `l` 像に入ること** | ☆**本ファイル**（進捗枠 **3 / 3** の部品が揃った、残りは配管 1） |

★位相の側（像が閉部分群）は第 772、群論の側（`α` と安定直線から `SL₂`）は
`§9-992` と `Lemma 3.1, (iv)` で済んでいる（`imageContainsSL2J_of_alpha'`）。

## ★★★★★★★★在庫（すでに `Found` にあるもの）

| 定理 | 内容 | 第 |
|---|---|---|
| `exists_sigma_smul_root_of_val` | `l ∤ v(q)` なら `σ(π) = η·π` なる `σ` が **`AdjoinRoot (Xˡ − C q)` の上に**ある | 994 |
| `not_lth_power_of_val` | `l ∤ v(q)` なら `q` は `l` 乗でない | 994 |
| `not_lth_power_of_not_dvd_val` | `l ∤ v_K(q)` なら `q` は `K(ζ_l)` でも `l` 乗でない（分岐指数は `l` と素） | 994 |
| `sigma_acts_as_alpha` | `σ(ζ)=ζ`・`σ(π)=ζπ` なら `σ(ζᵃπᵇ) = ζ^{a+b}πᵇ` | 993 |
| `upperM_one_mulVec` | `α = (1 1 / 0 1)` は `(a,b) ↦ (a+b, b)` | 993 |
| `tatePhi_pointMap` | Tate 一意化は `Point.map` と可換（同変性の群の形） | — |
| `tate_stab_of_pointStab` | 点の側の安定性を単数の側へ移す段 | — |

☆すなわち**両端は揃っている**——`Lˣ/qℤ` の側の `σ` と、行列 `α` の作用と、
Tate 一意化の同変性である。★足りないのは**それらを繋ぐ 3 本**である。

## ★★★★★★★★★★★★節点（進捗枠 **3 / 3** の部品が揃った）

| # | 節点 | 内容 | 重み |
|---|---|---|---|
| 1 | `tateTorsion_basis_zeta_pi` | ★**核はちょうど `lℤ × lℤ`（第 1161-1162、`zeta_pi_mem_zpowers_iff`）**。☆残るのは全射の側（`E[l]` が実際にこの 2 元で生成されること）と `Q` の無限位数の紐付け | ★★**閉じた**（第 1165 `zeta_pi_basis`） | 8 → **0** |
| 2 | `alpha_mem_localImage` | ★**座標の側（第 1164）と群論の側（第 1170 `exists_zeta_pi_of_torsion`）は揃った**。★★**第 1171 で完成**（`torsion_eq_phi_zeta_pi`） | 8 → **0** |
| 3 | `alpha_mem_globalImage` | ★**制限準同型は第 1167 で取れた**（`restrictLocalHom`・`restrictLocalHom_commutes`）。☆残るのは `L̄ ↪ M` を実際に取る段（`IsAlgClosed.lift`）と `galRep` との合成 | ★★**閉じた**（第 1168 `globalOfLocalHom`）。★★**第 1169 で完備**（`mem_map_of_mem_map_comp`） | 6 → **0** |

☆総重み 22 → **1**（第 1161-1165 で節点 1、第 1167-1169 で節点 3、
第 1170-1171 で節点 2 が閉じた）。

### ★★★★★★★★★★★★第 1171 で 3 節点がすべて揃った

`torsion_eq_phi_zeta_pi`——`Φ : Additive (G ⧸ ⟨Q⟩) ≃+ P` を Tate 一意化と読めば、
`l • Φ(c) = 0` なら `cˡ = 1`、したがって `c = [ζᵃπᵇ]`。
☆すなわち **`E[l]` の元はすべて `Φ [ζᵃπᵇ]` の形**である。

★★★**残るのは具体の `TateSetup`・`galRep` に当てはめる配管だけ**である。

### ★★★★★★★★第 1200——Kummer 拡大を Tate 設定の体として使える

`irreducible_X_pow_sub_C_of_not_pow`・`fact_irreducible_X_pow_sub_C_of_not_pow`
（`Found/GenEll/Thm38KummerExists.lean`）。

☆`l` 素数で `q` が `l` 乗でなければ `Xˡ − C q` は既約なので
`AdjoinRoot (Xˡ − C q)` は**体**になる。
★`TateSetup` は `[Field K]` を要求するので、これがないと当てはめられなかった。
☆`q` が `l` 乗でないことは `l ∤ v(q)` から出る（`not_lth_power_of_val`、第 994）。

### ★★★★★★★★★★★★★★★★第 1211——**単数の側は済んだ**

`exists_units_sigma_kummer`（`Found/GenEll/Thm38KummerExists.lean`）:

    l ∤ v(q) → ∃ π : Kˣ, ∃ σ : Kˣ →* Kˣ,
      πˡ = q ∧ σ ζ = ζ ∧ σ π = ζ · π     （K = AdjoinRoot (Xˡ − C q)）

☆これは `tate_sigma_coord_alpha`（第 1174）が**受け取る形そのもの**である
——第 994 の代数同型を `Units.map` で単数に落とし、
第 1178 の根 `π` と組み合わせた。
★`σ ζ = ζ` は `σ` が `K₀`-代数同型だから（`AlgEquiv.commutes`）。

★★★**残るのは `TateSetup` を `K` へ底変換する段だけ**である。

### ★★★★★★★★★★★★第 1172-1174——仮説がすべて消えた

| 仮説 | どこから出るか | 第 |
|---|---|---|
| `hQinf`（`Q` は無限位数） | `TateSetup` の `0 < v(q)` | 1172 |
| `hμ`（`μ_l = ⟨ζ⟩`） | `IsPrimitiveRoot`（整域） | 1173 |
| `hζl`・`hζprim` | `IsPrimitiveRoot`（単数群への持ち上げ） | 1174 |

★`tate_sigma_coord_alpha` と `tate_torsion_eq_phi_zeta_pi`（第 1174）が
**Tate 設定での完成形**であり、受けているのは
「`ζ` が原始 `l` 乗根」「`πˡ = q`」「`σ(ζ) = ζ`・`σ(π) = ζπ`」だけである。
☆残るのはこれらを `L_v(ζ_l, q^{1/l})` で実際に取る段（第 994 の Kummer の `σ`）と、
`galRep` の行列に読み替える段である。

### ★★★★★★★★第 1170 で節点 2 の群論の側が閉じた

`exists_zeta_pi_of_torsion`——`G ⧸ ⟨Q⟩` の `l`-捻れはすべて `[ζᵃπᵇ]` である。
☆Tate 一意化 `Φ : G ⧸ ⟨Q⟩ ≃+ E.Point` は加法同型なので、
`E[l]` の元は `Φ [ζᵃπᵇ]` の形に書ける。
★★**残るのは `Φ` を介して `E[l]` の言葉に直す配管だけ**である。
★★**残るのは具体の `TateSetup`・`galRep` への当てはめだけ**である。

### ★★★★★★★★第 1169 で節点 3 が完備した

`mem_map_of_mem_map_comp`——`f` を `galRep`、`g` を `globalOfLocalHom`、
`r` を `glRedPadic l` と読むと、これが
`imageContainsSL2J_of_alpha'` の `halpha` を**局所から受け取る形**である。

### ★★★★★★★★第 1168 で節点 3 が閉じた

`globalOfLocalHom`——`M` が代数閉体なら `L̄ = AlgebraicClosure L` は
`IsAlgClosed.lift` で `M` に埋め込めるので、その埋め込みで
`restrictLocalHom` を当てればよい。
★★`Gal(M/L_v) →* Gal(L̄/L)` ができた——**局所で実現される行列は大域の像にも入る**。

### ★★★★★★★★第 1167 で取れたもの（節点 3 の核）

`Found/GenEll/Thm38Decomposition.lean` の **`restrictLocalHom`**——
塔 `L → L_v → M` と `L` 上正規な中間体 `E` に対し

    `σ : M ≃ₐ[L_v] M`  →  `σ|_L : M ≃ₐ[L] M`  →  `σ|_E : E ≃ₐ[L] E`

の準同型。☆`restrictLocalHom_commutes` で埋め込みと可換であることも取った。

★★**単射性は要らない**——必要なのは「局所の元の像が大域の像に**含まれる**」
ことだけで、それは準同型があれば出る。

### ★★★★★★★★第 1161 で取れたもの

`Found/GenEll/Thm38ZetaPi.lean` の **`zeta_pi_indep`**——
`ζᵃ·πᵇ ∈ ⟨Q⟩` なら `l ∣ a` かつ `l ∣ b` である（★無条件）。

☆機構: `ζᵃπᵇ = Qⁿ` を `l` 乗すると `Qᵇ = Q^{ln}`、`Q` は無限位数なので `b = ln`。
すると `πᵇ = (πˡ)ⁿ = Qⁿ` なので `ζᵃ = 1`、`ζ` の位数は `l` なので `l ∣ a`。
★`orderOf_eq_of_primitive`（原始 `l` 乗根の位数はちょうど `l`）も一緒に取った。

★第 1162 で逆向き（`l ∣ a` かつ `l ∣ b` なら `ζᵃπᵇ ∈ ⟨Q⟩`）も取り、
**核がちょうど `lℤ × lℤ`** であることが iff の形で出た（`zeta_pi_mem_zpowers_iff`）。
☆これが「`([ζ], [π])` は `E[l]` の `ℤ/l`-基底である」の正確な形であり、
節点 2 が消費する界面である。

★第 1163 でさらに `zeta_pi_coord_eq_iff`——
`[ζ^{a₁}π^{b₁}] = [ζ^{a₂}π^{b₂}] ⟺ a₁ ≡ a₂, b₁ ≡ b₂ (mod l)`——を取った。
☆これで座標は `(ZMod l)²` の元として**一意に**決まる。
★`sigma_acts_as_alpha`（第 993）の `(a,b) ↦ (a+b, b)` と合わせれば、
`σ` の行列がちょうど `α = (1 1 / 0 1)` になる。

★★第 1164 でそれを組んだ `sigma_coord_alpha` を取った——
`σ(ζᵃπᵇ)` を座標 `(a', b')` で書いたなら必ず `a' ≡ a + b`、`b' ≡ b`（mod `l`）。
☆**これで座標の側は全部揃った**。

★★★第 1165 で**全射の側**（`zeta_pi_span`: `xˡ ∈ ⟨Q⟩` なら `x = ζᵃπᵇ`）も取り、
`zeta_pi_basis`（生成＋一意性）で**節点 1 が閉じた**。
☆残るのは Tate 一意化で `E[l]` を `Lˣ/qℤ` の `l`-捻れと同定する段（節点 2 の残り）と、
分解群経由で大域へ移す段（節点 3）だけである。

★節点 1-3 が閉じれば `imageContainsSL2J_of_alpha'` の `halpha` が埋まる。

### ★★★★節点 3 の道（測定）

`Gal(L̄_v/L_v) → Gal(L̄/L)` は**単射ではない**が、`L̄ ↪ L̄_v` を取れば
制限写像 `Gal(L̄_v/L_v) → Gal(L̄/L)` が定義され、その像が分解群である。
★`galRep` は `Gal(L̄/L)` の表現なので、局所の `σ` の像は大域の像に含まれる。
☆mathlib の `AlgHom.restrictNormal` / `IsAlgClosed.lift` が素材である。

### ☆節点 1 が本体である理由

`E[l]` の `ℤ/l`-基底を `(ζ, π)` に取ることが、`α = (1 1 / 0 1)` という**具体的な行列**を
出す唯一の道である。★`Lˣ/qℤ` の `l`-捩れは `μ_l · π^ℤ / qℤ` であり、
`πˡ = q` なので `π` の類は位数 `l`、`ζ` の類も位数 `l`、両者は独立である。

## ★★★★★★何が `Theorem 3.8` に残るか（第 1160 の総括）

    I（`l`-巡回の読み、残り 4） ＋ II（本ファイル、残り 1） ＝ **5**

☆`Theorem 3.8` 以外（`Corollary 4.3`・`4.4`）は本定理の上にしか立たない。
★`Theorem 2.1`（§2）は曲線の Riemann–Roch が mathlib に無いので別筋である。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta

/-! ## ★出典の紐付け(`.src`) -/

def alphaBridgeFrame.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(α が mod l 像に入る段の橋——3 節点の枠)",
    sectionId := "genell-thm-3-8" }

def alphaBridgeFrame.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_sigma_smul_root_of_val(σ の存在、第 994、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_sigma_smul_root_of_val") 1,
    .citation "[ABC3]" "sigma_acts_as_alpha(行列表示、第 993、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.sigma_acts_as_alpha") 1,
    .citation "[ABC3]" "imageContainsSL2J_of_alpha'(残る仮説は α だけ、第 772、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.imageContainsSL2J_of_alpha'") 1,
    .implicitStep
      ("★★★★**2026-09-01（第 1160）の測定**——`Theorem 3.8` を `Found` で閉じるのに" ++
       "要るのはちょうど 2 つ: (I) `l`-巡回の読みの橋（`LCyclicReading`、残り重み 5）と " ++
       "(II) `α` が mod `l` 像に入ること（本ファイル、第 1161 で 22 → 17）。" ++
       "☆両端——`Lˣ/qℤ` の側の `σ`（第 994）と行列 `α` の作用（第 993）と " ++
       "Tate 一意化の同変性（`tatePhi_pointMap`）——はすべて揃っている。" ++
       "★足りないのは (1) `E[l]` の基底を `(ζ, π)` に取る段、" ++
       "(2) 局所の mod `l` 像が `α` を含むこと、" ++
       "(3) 分解群経由で大域の像へ移す段の 3 本である。") 22 ]

end ABC3.Skeleton.GenEll
